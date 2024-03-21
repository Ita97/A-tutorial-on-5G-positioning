clc;  close all;
clearvars -except viewer;

%% Simulation Variables
rng(42);                       % Set RNG state for repeatability
pd = truncate(makedist("Normal", "sigma",50),-100,100); % probability distribution
                            
reflections_order = 2;       % 0 (LOS) default | 1 (NLOS)
plot_rays = true;            % plot raytracer rays on map
refinement = true;          % true (P2) default | false (only P1)

num_tests = 1;             
dimensions = 3;            % 2/3-D
K = 1000;                  % number of iteration for NLS algorithm
stop_cond = 1e-4;

drop_nlos = true;           % drop NLOS BS
BS_synch = true;            % BS synchronized?
max_gNBs = 5;               % Max number of BS used

%% Simulation Parameters
% Bandwidth configuration, required to set the channel sampling rate and for perfect channel estimation
% (52 RBs at 15 kHz SCS for 10 MHz BW)
% (106 RBs at 15 kHz SCS for 20 MHz BW)
% (216 RBs at 15 kHz SCS for 40 MHz BW)
% (270 RBs at 15 kHz SCS for 50 MHz BW)

% (51 RBs at 30 kHz SCS for 20 MHz BW)
% (78 RBs at 30 kHz SCS for 30 MHz BW)
% (106 RBs at 30 kHz SCS for 40 MHz BW)
% (217 RBs at 30 kHz SCS for 80 MHz BW)
% (273 RBs at 30 kHz SCS for 100 MHz BW)

% (24 RBs at 60 kHz SCS for 20 MHz BW)
% (51 RBs at 60 kHz SCS for 40 MHz BW)
% (135 RBs at 60 kHz SCS for 100 MHz BW)
% (273 RBs at 60 kHz SCS for 200 MHz BW)

% physical layer
nFrames = 1;
SCS = 30;                             % sub-carrier spacing (kHz)
NSizeGrid = 273;                      % number of RBs

NSubCar= 12;
fc = 3.8e9;                           % carrier frequency (Hz)
BW = SCS*NSizeGrid*NSubCar*1000;           % BandWidth (Hz)

if SCS<60 || ( SCS==60 && NSizeGrid<=135)
    frequencyRange = 'FR1';
    NoiseFigureRx = 9;          % Noise Figure @ UE [dB]
    NoiseFigureTx = 5;          % Noise Figure @ BS [dB]
    %fc = 3.8e9;
    if ~(fc < 7.125e9)
        error('In FR1 central frequency must be below 7.125GHz')
    end
    disp('Frequency Range 1');
else
    frequencyRange = 'FR2';
    fc = 27.7e9;
    NoiseFigureRx = 10;          % Noise Figure @ UE [dB]
    NoiseFigureTx = 7;           % Noise Figure @ BS [dB]
    if ~(fc > 24.250e9 && fc < 52.6e9)
        error('In FR2 central frequency must be between 24.250GHz and 52.600GHz')
    end
    disp('Frequency Range 2');
end
disp(['Bandwidth ' num2str(BW/1000000) ' MHz'])

% physical constant 
c = physconst('LightSpeed');          % speed of light
lambda = c/fc;                        % wave length

% Rx/Tx antennas values
T_ant = 298;                % Antenna temperature [K] (25°C)
TxPower = 33;               % BS Tx Power [dBm] Urban
%TxPower = 49;               % BS Tx Power [dBm] Rural
RxPower = 23;               % UE Tx Power [dBm]

ref_site = 1;               % BS used as reference for geo<->cartesian conversion
TxArraySize = [4 4];        % BS Antenna Array Size
RxArraySize = [2 2];        % UE Antenna Array Size 
numTxAnt = prod(TxArraySize);
numRxAnt = prod(RxArraySize);
TxAzimuthRange = [-60 60];
TxElevationRange = [-90 0];
RxAzimuthRange = [-180 180];
RxElevationRange = [0 90];
ueArrayOrientation = [0 45].';      % azimuth (0 deg is E, 90 deg is N) and 
                                    % elevation (positive points upwards) in deg


%% UE positions
load ue_positions/giuriati.mat
numUEPos = size(listUEPos,1); 

%% gNBs positions
n_cell = 3; % number of cells
h_gNB = 4; % height [m]
load bs_positions/campus_leonardo.mat % bs positions and map name

numBSs = size(gNBPos,1);
TxArrayOrientation = [[0 0]; [120 0]; [-120 0]];
gNBs = cell(numBSs, 1);
pci = [];
for gNBIdx = 1:numBSs
    for i = 1:n_cell
        t_pci = randperm(1008,1) - 1;
        while ismember(t_pci,pci)
            t_pci = randperm(1008,1) - 1;
        end
        cells(i) = gNBCellClass(gNBPos(gNBIdx,:),t_pci,TxArrayOrientation(i,:), TxArraySize, fc);
        pci((gNBIdx-1)*n_cell+i) = t_pci;
    end
    gNBs{gNBIdx} = cells;
end
numgNBs = size(gNBs,1);

%% Setup Viewer

if exist('viewer','var') && isvalid(viewer) % viewer handle exists and viewer window is open
    viewer.clearMap();
else
    viewer = siteviewer("Basemap","openstreetmap","Buildings",map);   
end

bsSite = cell(1,numgNBs);
bsArray = cell(1,numgNBs);

for gNBIdx = 1 : numgNBs 
    bsSite{gNBIdx} = txsite("Name","BS"+gNBIdx, ...
        "Latitude",gNBs{gNBIdx}(:).getCoordinates(1),"Longitude",gNBs{gNBIdx}(:).getCoordinates(2),...
        "AntennaAngle",gNBs{gNBIdx}(:).getOrientation()',...
        "AntennaHeight",h_gNB,...  % in m
        "TransmitterFrequency",fc, ...
        "TransmitterPower",db2pow(TxPower-30),...
        "Antenna",{gNBs{gNBIdx}(:).antennaArray}); % in W
    show(bsSite{gNBIdx})
end

ueSite = cell(1,numUEPos);
ueArray = cell(1,numUEPos);
for UEIdx = 1 : numUEPos
    ueSite{UEIdx} = rxsite("Name","UE"+UEIdx, ...
            "Latitude",listUEPos(UEIdx,1),"Longitude",listUEPos(UEIdx,2),...
            "AntennaHeight",2,... % in m
            "AntennaAngle",ueArrayOrientation(1:2));
%     ueArray{UEIdx} = phased.URA('Size',[RxArraySize(1:2) 1 1],'Spacing', [0.5*lambda*[1 1] 1 1]);
    ueArray{UEIdx} = phased.NRRectangularPanelArray('Size',[RxArraySize(1:2) 1 1],'Spacing', [0.5*lambda*[1 1] 1 1]);
    ueArray{UEIdx}.ElementSet = {phased.IsotropicAntennaElement};   % isotropic antenna element
end
   
%% Convert geographic coordinates to cartesian

list_UE_xyz = zeros(numUEPos, 3);
for UEIdx = 1:numUEPos
    list_UE_xyz(UEIdx,:) = from_lat_long_to_xyz(ueSite{UEIdx}, ueSite{ref_site}); % position to convert and position reference 
end

for gNBIdx = 1:numgNBs
    mat_xyz =from_lat_long_to_xyz(bsSite{gNBIdx}, ueSite{ref_site});% position to convert and position reference
%     gNBs{gNBIdx}(:).xyz_cartesian = mat_xyz;
    for i = 1:n_cell
        gNBs{gNBIdx}(i).xyz_cartesian = mat_xyz(i,:);
    end
end

%% gNB Carrier
% Configure carrier properties
carrier = repmat(nrCarrierConfig("NSizeGrid", NSizeGrid,"SubcarrierSpacing", SCS),numgNBs,n_cell);
for gNBIdx = 1:numgNBs
    
    for cellIdx = 1:n_cell
        carrier(gNBIdx,cellIdx).NCellID = gNBs{gNBIdx}(cellIdx).pci;
    end
end
validateCarriers(carrier);

%% UE Carrier
ue = nrCarrierConfig;
% Configure carrier properties
ue.NSizeGrid = NSizeGrid;
ue.SubcarrierSpacing = SCS;
validateCarriers(ue);
ofdmInfo = nrOFDMInfo(ue);

%% Noise

NF = 10^(NoiseFigureRx/10);  % Noise Figure 
Te = T_ant + 290*(NF-1);     % noise temperature [K]
kB = physconst("Boltzmann"); % [J/K]
N0DL = sqrt((kB*ofdmInfo.SampleRate*Te)/2);     % noise amplitude DL

% Shadowing Fading (SF)
sigma_SF = 3; % dB

%% Beamforming Plotting parameters
prm.TxArraySize = TxArraySize;       % Transmit array size, [rows cols]
prm.TxAZlim = [-60 60];             % Transmit azimuthal sweep limits
prm.TxELlim = [-90 0];              % Transmit elevation sweep limits

prm.RxArraySize = RxArraySize;        % Receive array size, [rows cols]
prm.RxAZlim = [-180 180];           % Receive azimuthal sweep limits
prm.RxELlim = [0 90];               % Receive elevation sweep limits

prm.ElevationSweep = frequencyRange=="FR2";         % Enable/disable elevation sweep
prm.RSRPMode = 'SSSwDMRS';          % {'SSSwDMRS', 'SSSonly'}
prm.NCellID = 1;                    % Cell ID
prm.FreqRange = frequencyRange;     % Frequency range: 'FR1' or 'FR2'
prm.CenterFreq = fc;                % Hz
ssbpattern = 'C';
if fc<=6e9
    switch SCS
        case 15
            ssbpattern = 'A';
        case 30
            ssbpattern = 'B'; % B or C
    end
else
    switch SCS
        case 120
            ssbpattern = 'D';
        case 240
            ssbpattern = 'E';
    end
end

prm.SSBlockPattern = ['Case ' ssbpattern];      % Case A/B/C/D/E
if fc < 3e9 
    if (ssbpattern=='A' || ssbpattern=='B' || ssbpattern=='C')
        ssbsymb=4;
    else
        ssbsymb=8;
    end
else
    if (ssbpattern=='A' || ssbpattern=='B' || ssbpattern=='C')
        ssbsymb=8;
    else
        ssbsymb=64;
    end
end

prm.SSBTransmitted = [ones(1,ssbsymb) zeros(1,0)];   % 4/8 or 64 in length
prm = validateParams(prm);
        
%% For each UE position
for posIdx = 1: numUEPos
    numgNBs = size(gNBs,1);
    UEPos = listUEPos(posIdx,:);
    UE_xyz = list_UE_xyz(posIdx, :); 
    cannot_reach = zeros(numgNBs, n_cell);
    reachable_gNBs_idx = [];
    reachable_gNBs =[];

    show(ueSite{posIdx});

    %% Ray Tracing Analysis
    % Perform ray tracing analysis using the Shooting and Bouncing Rays (SBR) method. 
    % The SBR method includes effects from surface reflections and diffractions but 
    % does not include effects from refraction or scattering.
    pm = propagationModel("raytracing", "Method", "sbr","MaxNumReflections",reflections_order);
    
    channel = cell(numgNBs,n_cell);
    scatPos=cell(numgNBs,n_cell);
    LOSRay = cell(numgNBs,n_cell);
    los = zeros(1, numgNBs);
    for gNBIdx = 1:numgNBs
            rays = raytrace(bsSite{gNBIdx},ueSite{posIdx},pm,"Type","pathloss");
            for cellIdx = 1:n_cell 
                % if the there are no rays connecting the gNBCell to the UE, skip it
                if isempty(rays{cellIdx})
                    cannot_reach(gNBIdx, cellIdx) = 1;
                    continue
                end
                num_rays = length([rays{cellIdx}.PropagationDelay]);
                %azlim = gNBs{gNBIdx}(cellIdx).getAzimuthLimits;
                ray_index = false(1,num_rays);
                for i = 1:num_rays
                    aod = rays{cellIdx}(i).AngleOfDeparture;
                    if gNBs{gNBIdx}(cellIdx).isInAzimuthRange(aod(1))
                        ray_index(i) = true;
                    end
                end
                if ~sum(ray_index)  %if no ray is belonging to the cell, skip it
                    cannot_reach(gNBIdx, cellIdx) = 1;
                    continue;
                end

                reachable_gNBs_idx = [reachable_gNBs_idx; [gNBIdx,cellIdx]];
                reachable_gNBs = [reachable_gNBs; gNBs{gNBIdx}(cellIdx)];
                cellRays = rays{cellIdx}(ray_index);
                
                N = size(cellRays,2);
                scatPos{gNBIdx,cellIdx} = [];
                SPSite = cell(1,N);
                for i = 1: N
                    M = cellRays(i).NumInteractions;
                    scatter = zeros(2,M);
                    for j = 1:M
                        scatter(:,j) = cellRays(i).Interactions(j).Location(1:2)';
                    end
                    if M>0
                        SPSite{i} = rxsite("Name","Scatter"+i, ...
                        "Latitude",scatter(1,:),"Longitude",scatter(2,:));
                    else
                        SPSite{i} = rxsite("Name","Scatter"+i, ...
                        "Latitude",0,"Longitude",0);
                    end
                    scatterPos_xyz = from_lat_long_to_xyz(SPSite{i}, ueSite{ref_site});
                    scatPos{gNBIdx,cellIdx} = [scatPos{gNBIdx,cellIdx};scatterPos_xyz];
                
                    LOSRay{gNBIdx,cellIdx}(i) = cellRays(i).LineOfSight;
                end
                
                pathToAs = [cellRays.PropagationDelay];                                  % Time of arrival of each ray
                avgPathGains  = -[cellRays.PathLoss];                                    % Average path gains of each ray
                pathAoDs = [cellRays.AngleOfDeparture];                                  % AoD of each ray
                pathAoAs = [cellRays.AngleOfArrival];                                    % AoA of each ray
                isLOS = any([cellRays.LineOfSight]);                                     % Line of sight flag

                if ~los(gNBIdx) && isLOS
                    los(gNBIdx) = isLOS;
                end
            
                if ~isempty(cellRays) && plot_rays
                    plot(cellRays)
                end 
                
                % Set Up CDL Channel Model
                channel{gNBIdx, cellIdx} = nrCDLChannel;
                channel{gNBIdx, cellIdx}.DelayProfile = 'Custom';
                channel{gNBIdx, cellIdx}.PathDelays = pathToAs;           % Time of arrival of each ray (normalized to 0 sec)
                channel{gNBIdx, cellIdx}.AveragePathGains = avgPathGains;
                channel{gNBIdx, cellIdx}.AnglesAoD = pathAoDs(1,:);       % azimuth of departure
                channel{gNBIdx, cellIdx}.AnglesZoD = 90-pathAoDs(2,:);    % channel uses zenith angle, rays use elevation
                channel{gNBIdx, cellIdx}.AnglesAoA = pathAoAs(1,:);       % azimuth of arrival
                channel{gNBIdx, cellIdx}.AnglesZoA = 90-pathAoAs(2,:);    % channel uses zenith angle, rays use elevation
                channel{gNBIdx, cellIdx}.HasLOSCluster = isLOS;           % Line of sight flag
                channel{gNBIdx, cellIdx}.CarrierFrequency = fc;
                channel{gNBIdx, cellIdx}.NormalizeChannelOutputs = false; % do not normalize by the number of receive antennas, this would change the receive power
                channel{gNBIdx, cellIdx}.NormalizePathGains = false;      % set to false to retain the path gains
                
                channel{gNBIdx, cellIdx}.ReceiveAntennaArray = ueArray{posIdx};
                channel{gNBIdx, cellIdx}.ReceiveArrayOrientation = [ueArrayOrientation(1); (-1)*ueArrayOrientation(2); 0];  % the (-1) converts elevation to downtilt
                
                channel{gNBIdx, cellIdx}.TransmitAntennaArray = gNBs{gNBIdx}(cellIdx).antennaArray;
                channel{gNBIdx, cellIdx}.TransmitArrayOrientation = [gNBs{gNBIdx}(cellIdx).orientation(1); (-1)*gNBs{gNBIdx}(cellIdx).orientation(2); 0];   % the (-1) converts elevation to downtilt 
            
                channel{gNBIdx, cellIdx}.SampleRate = ofdmInfo.SampleRate;
            end
        
    
    end
        
    if numel(logical(los)) ~= numgNBs
            error('nr5g:InvalidLOSLength',['Length of line of sight flag (' num2str(numel(logical(los))) ...
            ') must be equal to the number of configured gNBs (' num2str(numgNBs) ').']);
    end
    numgNBs = numel(reachable_gNBs);        % we consider from now on, only the reachable cell

    %% Drop NLOS gNBs
    if drop_nlos
        LOS = false(numgNBs, 1);
        for i= 1: numgNBs
            LOS(i) = channel{reachable_gNBs_idx(i,1), reachable_gNBs_idx(i,2)}.HasLOSCluster;
        end
        if any(not(LOS))
            reachable_gNBs_idx = reachable_gNBs_idx(LOS,:);
            reachable_gNBs = reachable_gNBs(LOS);
            numgNBs = numel(reachable_gNBs);
        end
    end

    disp([13 'For position ' num2str(posIdx) ' the number of available gNB cells are: ' num2str(numgNBs)])

    cell_carrier = cell(numgNBs,1);
    for i = 1:numgNBs
            cell_carrier{i} = carrier(reachable_gNBs_idx(i,1), reachable_gNBs_idx(i,2));
    end

    

    %% PRS definition 
    % for P-2
    prsSym = cell(1,numgNBs);
    prsInd = cell(1,numgNBs);

    % for each gNBs assign a different Slot
    gap = 2;
    prsSlotOffsets = 0:gap:(gap*numgNBs - 1);
    prsIDs = randperm(4096,numgNBs) - 1;
    
    % Configure PRS properties
    numPRSRes = 12; % one for each beam
    prs = nrPRSConfig;
    prs.PRSResourceSetPeriod = [12 0];
    prs.PRSResourceOffset = ones(1,numPRSRes);
    prs.REOffset = 0:(numPRSRes-1);
    prs.PRSResourceRepetition = 1;
    prs.PRSResourceTimeGap = 1;
    prs.MutingPattern1 = [];
    prs.MutingPattern2 = [];
    prs.NumRB = NSizeGrid;
    % numRB = prs.NumRB;
    prs.RBOffset = 0;
    prs.CombSize = 12;
    prs.NumPRSSymbols = 12*ones(1, numPRSRes);
    % numPRSSymbol=prs.NumPRSSymbols;
    prs.SymbolStart = zeros(1,numPRSRes);
    prs = repmat(prs,1,numgNBs);
    for gNBIdx = 1:numgNBs
        prs(gNBIdx).PRSResourceOffset = prsSlotOffsets(gNBIdx)*ones(1,numPRSRes);
        prs(gNBIdx).NPRSID = prsIDs(gNBIdx);
        cell_carrier{gNBIdx}.NSlot = prs(gNBIdx).PRSResourceSetPeriod(2)+prs(gNBIdx).PRSResourceOffset(1);
        tmpSym = nrPRS(cell_carrier{gNBIdx},prs(gNBIdx));
        tmpInd = nrPRSIndices(cell_carrier{gNBIdx},prs(gNBIdx));
        prsSym{gNBIdx}=nrPRS(cell_carrier{gNBIdx},prs(gNBIdx), "OutputResourceFormat","cell");
        prsInd{gNBIdx}=nrPRSIndices(cell_carrier{gNBIdx},prs(gNBIdx), "OutputResourceFormat","cell");
    end
    powerPRS = 0;
    ports=1;

    
    %% Synchronization Signal Burst Configuration and Generation

    txBurst = nrWavegenSSBurstConfig;
    txBurst.BlockPattern = prm.SSBlockPattern;
    txBurst.TransmittedBlocks = prm.SSBTransmitted;
    txBurst.Period = 20;
%             txBurst.NCRBSSB = 10;
%             txBurst.KSSB = 2*4; % for D even and for E multiple of 4
    txBurst.SubcarrierSpacingCommon = prm.SubcarrierSpacingCommon;
    
    % Configure an nrDLCarrierConfig object to use the synchronization signal
    % burst parameters and to disable other channels. This object will be used
    % by nrWaveformGenerator to generate the SS burst waveform.
    
    cfgDL = configureWaveformGenerator(prm,txBurst);
%             cfgDL.SampleRate=ofdmInfo.SampleRate;
    burstWaveform = nrWaveformGenerator(cfgDL); 
    
    % SSB Carrier
    ssb_carrier = cell(1,numgNBs);
    for i = 1:numgNBs
        ssb_carrier{i} = cell_carrier{i};
        ssb_carrier{i}.NSizeGrid = cfgDL.SCSCarriers{1}.NSizeGrid;
    end
    ssb_ofdmInfo = nrOFDMInfo(ssb_carrier{1});

    %% For each Test
    for test = 1:num_tests

        disp(['Test #' num2str(test)]);
        numgNBs = length(reachable_gNBs);
        u_0 = randi([-100,100],1,dimensions);
        if dimensions==3
            u_0(1,3) = 0;
        end
        
%         if test==1
%             % Scatter Plot
%             tmpScatter = cell(1,numgNBs);
%             tmpRayLOS = cell(1,numgNBs);
%             for i = 1:numgNBs
%                 tmpScatter{i} = scatPos{reachable_gNBs_idx(i,1), reachable_gNBs_idx(i,2)};
%                 tmpRayLOS{i} = LOSRay{reachable_gNBs_idx(i,1), reachable_gNBs_idx(i,2)};
%             end
%             scatPos = tmpScatter;
%             LOSRay = tmpRayLOS;
%             clear tmpRayLOS;
%             clear tmpScatter;
%         end
        %% Add Signal Delays
    
        sampleDelay = cell(1,numgNBs);
        cell_channel = cell(1,numgNBs);
        for i = 1:numgNBs
            cell_channel{i} = channel{reachable_gNBs_idx(i,1),reachable_gNBs_idx(i,2)};
            chInfo = info(cell_channel{i});
            sampleDelay{i} = ceil((chInfo.PathDelays)*cell_channel{i}.SampleRate);
        end
        maxChDelay = max(cell2mat(sampleDelay))+ chInfo.ChannelFilterDelay;

        if test>1 % return back the channel order
            cell_carrier = cell(numgNBs,1);
            cell_channel = cell(numgNBs,1);
            for i = 1:numgNBs
                idx = [ reachable_gNBs_idx(i,1), reachable_gNBs_idx(i,2) ];
                cell_carrier{i} = carrier(idx(1), idx(2));
                cell_channel{i} = channel{idx(1), idx(2)};
                reset(cell_channel{i});
            end
        end

        %% Transmit-End Beam Sweeping
        % To achieve TRP beam sweeping, beamform each of the SS blocks in the generated burst using analog beamforming. 
        % Based on the number of SS blocks in the burst and the sweep ranges specified, 
        % determine both the azimuth and elevation directions for the different beams. 
        % Then beamform the individual blocks within the burst to each of these directions.
        
        % Number of beams at both transmit and receive ends
        numBeams = sum(txBurst.TransmittedBlocks);
        
        % Transmit beam angles in azimuth and elevation, equi-spaced
        txBeamAng = cell(1, numgNBs);

        for gNBIdx = 1 : numgNBs
            azTxBW = beamwidth(reachable_gNBs(gNBIdx).antennaArray,fc,'Cut','Azimuth');
            elTxBW = beamwidth(reachable_gNBs(gNBIdx).antennaArray,fc,'Cut','Elevation');
            txBeamAng{gNBIdx} = hGetBeamSweepAngles(numBeams,prm.TxAZlim,prm.TxELlim, ...
                azTxBW,elTxBW,prm.ElevationSweep);
            % For evaluating transmit-side steering weights
            SteerVecTx = phased.SteeringVector('SensorArray',reachable_gNBs(gNBIdx).antennaArray, ...
                'PropagationSpeed',c);
        end
              
        % Get the set of OFDM symbols occupied by each SSB
        numBlocks = length(txBurst.TransmittedBlocks);
        burstStartSymbols = ssBurstStartSymbols(txBurst.BlockPattern,numBlocks);
        burstStartSymbols = burstStartSymbols(txBurst.TransmittedBlocks==1);
        burstOccupiedSymbols = burstStartSymbols.' + (1:4);
        
        % Apply steering per OFDM symbol for each SSB
        gridSymLengths = repmat(ssb_ofdmInfo.SymbolLengths,1,cfgDL.NumSubframes);
        %   repeat burst over numTx to prepare for steering
        strTxWaveform = cell(1, numgNBs);
        wT = cell(1, numgNBs);
        for gNBIdx = 1 : numgNBs
            strTxWaveform{gNBIdx} = repmat(burstWaveform,1,prm.NumTx)./sqrt(prm.NumTx);
        end
        
        for ssb = 1:numBeams
            
            % Extract SSB waveform from burst
            blockSymbols = burstOccupiedSymbols(ssb,:);
            startSSBInd = sum(gridSymLengths(1:blockSymbols(1)-1))+1;
            endSSBInd = sum(gridSymLengths(1:blockSymbols(4)));
            for gNBIdx = 1 : numgNBs
                ssbWaveform = strTxWaveform{gNBIdx}(startSSBInd:endSSBInd,1);
            
                % Generate weights for steered direction
                wT{gNBIdx} = SteerVecTx(fc,txBeamAng{gNBIdx}(:,ssb));
            
                % Apply weights per transmit element to SSB
                strTxWaveform{gNBIdx}(startSSBInd:endSSBInd,:) = ssbWaveform.*(wT{gNBIdx}'); 
            end
        end

        ibar_SSB = 0;
        pbchdmrsRef = nrPBCHDMRS(ssb_carrier{1}.NCellID,ibar_SSB);
        pbchDMRSInd = nrPBCHDMRSIndices(ssb_carrier{1}.NCellID);
        pssRef = nrPSS(ssb_carrier{1}.NCellID);
        pssInd = nrPSSIndices;
        pssGrid = zeros([240 4]);
        pssGrid(pssInd) = pssRef;
        pssGrid(pbchDMRSInd) = pbchdmrsRef;
        refGrid = zeros([12*ssb_carrier{1}.NSizeGrid ssb_ofdmInfo.SymbolsPerSlot]);
        burstOccupiedSubcarriers = ssb_carrier{1}.NSizeGrid*6 + (-119:120).';
        refGrid(burstOccupiedSubcarriers, burstOccupiedSymbols(1,:)) = pssGrid;

        %% Receive-End Beam Sweeping and Measurement
        % Receive beam angles in azimuth and elevation, equi-spaced
        azRxBW = beamwidth(ueArray{posIdx},fc,'Cut','Azimuth');
        elRxBW = beamwidth(ueArray{posIdx},fc,'Cut','Elevation');
        rxBeamAng = hGetBeamSweepAngles(numBeams,prm.RxAZlim,prm.RxELlim, ...
            azRxBW,elRxBW,prm.ElevationSweep);
        
        % For evaluating receive-side steering weights
        SteerVecRx = phased.SteeringVector('SensorArray',ueArray{posIdx}, ...
            'PropagationSpeed',c);

        %% Loop over all receive beams
        rsrp = cell(1, numgNBs);
        wR = cell(1, numgNBs);
        for gNBIdx = 1:numgNBs
            reset(cell_channel{gNBIdx});
            rsrp{gNBIdx} = zeros(numBeams,numBeams);
            spLoss = -(cell_channel{gNBIdx}.AveragePathGains);
            rxGain = 10.^(spLoss/20);
            % Fading channel, with path loss
            txWave = [strTxWaveform{gNBIdx}; zeros(maxChDelay,size(strTxWaveform{gNBIdx},2))];
            fadWave = cell_channel{gNBIdx}(txWave);
        
            % Receive gain, to compensate for the path loss
            fadWaveG = fadWave*min(rxGain);
        
            % Add WGN
            noise = N0DL*complex(randn(size(fadWaveG)),randn(size(fadWaveG)));
            
            rxWaveform = fadWaveG + noise;
%                 rxWaveform = fadWave;
            for rIdx = 1:numBeams  
            
                % Generate weights for steered direction
                wR{gNBIdx} = SteerVecRx(fc,rxBeamAng(:,rIdx));
            
                % Apply weights per receive element
                if strcmp(frequencyRange, 'FR1')
                    strRxWaveform = rxWaveform.*(wR{gNBIdx}');
                else  % for FR2, combine signal from antenna elements
                    strRxWaveform = rxWaveform*conj(wR{gNBIdx});
                end
            
                % Correct timing
                offset = nrTimingEstimate(ssb_carrier{gNBIdx}, ...
                    strRxWaveform(1:ssb_ofdmInfo.SampleRate*1e-3,:),refGrid*wR{gNBIdx}(1)');
                if offset > maxChDelay
                    offset = 0;
                end
                strRxWaveformS = strRxWaveform(1+offset:end,:);
            
                % OFDM Demodulate
                rxGrid = nrOFDMDemodulate(ssb_carrier{gNBIdx},strRxWaveformS);

            
                % Loop over all SSBs in rxGrid (transmit end)
                for tIdx = 1:numBeams
                    % Get each SSB grid
                    rxSSBGrid = rxGrid(burstOccupiedSubcarriers, ...
                        burstOccupiedSymbols(tIdx,:),:);
                    if strcmpi(prm.RSRPMode,'SSSwDMRS')
                        meas = nrSSBMeasurements(rxSSBGrid,ssb_carrier{gNBIdx}.NCellID,mod(tIdx-1,8));
                    else
                        meas = nrSSBMeasurements(rxSSBGrid,ssb_carrier{gNBIdx}.NCellID);
                    end
                    % Make measurements, store per receive, transmit beam
%                         rsrp{gNBIdx}(rIdx,tIdx) = measureSSB(rxSSBGrid,prm.RSRPMode,ssb_carrier{gNBIdx}.NCellID);
                    rsrp{gNBIdx}(rIdx,tIdx) = max(meas.RSRPPerAntenna);     
                
                end
            end
        end
        %% Beam Determination
        %  determine the best beam-pair link based on the RSRP measurement.
        coarse_AoD = zeros(numgNBs, 2);
        coarse_AoA = zeros(numgNBs, 2);
        refined_AoD = zeros(numgNBs, 2);
        refined_AoA = zeros(numgNBs, 2);
        s_aoa = [];

        for gNBIdx = 1: numgNBs
            s_aoa = [s_aoa; reachable_gNBs(gNBIdx).getXYZ()];
            [m,i] = max(rsrp{gNBIdx},[],'all','linear');    % First occurrence is output
            % i is column-down first (for receive), then across columns (for transmit)
            [rxBeamID,txBeamID] = ind2sub([numBeams numBeams],i(1));
            wT{gNBIdx} = SteerVecTx(fc,txBeamAng{gNBIdx}(:,txBeamID));
            wR{gNBIdx} = SteerVecRx(fc,rxBeamAng(:,rxBeamID));
        
            coarse_AoD(gNBIdx, :) = txBeamAng{gNBIdx}(:,txBeamID);
            coarse_AoA(gNBIdx, :) = rxBeamAng(:,rxBeamID);
%                 coarse_AoD(gNBIdx, :) = txBeamAng{gNBIdx}(:,txBeamID)+reachable_gNBs(gNBIdx,:).orientation';
%                 coarse_AoA(gNBIdx, :) = rxBeamAng(:,rxBeamID)+ueArrayOrientation;
            if abs(coarse_AoA(gNBIdx,1)) > 180
                coarse_AoA(gNBIdx,1) = coarse_AoA(gNBIdx,1)-sign(coarse_AoA(gNBIdx,1))*360;
            end
            if abs(coarse_AoA(gNBIdx,2)) > 90
                coarse_AoA(gNBIdx,2) = coarse_AoA(gNBIdx,2)-sign(coarse_AoA(gNBIdx,2))*180;
            end

            if plot_rays
                % Plot Radiation Patterns
                % Plot the radiation patterns obtained for the UE and the base station.
                nLayers = 1;
                scOffset = 0;   % no offset
                noRBs = cfgDL.SCSCarriers{1}.NSizeGrid;      % average channel conditions over 1 RB to calculate beamforming weights
                % Plot UE radiation pattern
                ueSite{posIdx}.Antenna = clone(cell_channel{gNBIdx}.ReceiveAntennaArray); % need a clone, otherwise setting the Taper weights would affect the channel array
                ueSite{posIdx}.Antenna.Taper = conj(wR{gNBIdx});
                pattern(ueSite{posIdx},fc,"Size",7);
                
                % Plot BS radiation pattern
                bsSite{reachable_gNBs_idx(gNBIdx,1)}(reachable_gNBs_idx(gNBIdx,2)).Antenna = clone(reachable_gNBs(gNBIdx).antennaArray); % need a clone, otherwise setting the Taper weights would affect the channel array
                bsSite{reachable_gNBs_idx(gNBIdx,1)}(reachable_gNBs_idx(gNBIdx,2)).Antenna.Taper = conj(wT{gNBIdx}); 
                pattern(bsSite{reachable_gNBs_idx(gNBIdx,1)}(reachable_gNBs_idx(gNBIdx,2)),fc,"Size",10);
            end

            % Get the azimuthal sweep range based on the SSB transmit beam direction
            % and its beamwidth in azimuth plane
            azSweepRange = [coarse_AoD(gNBIdx, 1) - azTxBW/2 coarse_AoD(gNBIdx, 1) + azTxBW/2];

            % Get the elevation sweep range based on the SSB transmit beam direction
            % and its beamwidth in elevation plane
            elSweepRange = [coarse_AoD(gNBIdx, 2) - elTxBW/2 coarse_AoD(gNBIdx, 2) + elTxBW/2];
            numPRSBeams = numPRSRes; %     TODO control if slot are muted
            prsTransmitted = ones(1, numPRSBeams);

            prsBeamAng = hGetBeamSweepAngles(numPRSBeams,azSweepRange,elSweepRange,azTxBW,elTxBW);


            if refinement
                % Calculate the steering vectors for all active PRS resources.
                wT{gNBIdx} = zeros(numTxAnt,numPRSBeams);
                
                for beamIdx = 1:numPRSBeams
                    wT{gNBIdx}(:,beamIdx) = SteerVecTx(fc,prsBeamAng(:,beamIdx));
                end 
                % Apply Digital Beamforming
                % Loop over all PRS resources and apply the digital beamforming to 
                % all the active ones. Digital beamforming is considered to offer frequency 
                % selective beamforming within the same OFDM symbol.
                
                % Initialize the beamformed grid
                bfGrid = nrResourceGrid(cell_carrier{gNBIdx},numTxAnt);
                % Get the active PRS resource indices
                activeRes = find(logical(prsTransmitted));
                for resIdx = 1:numPRSBeams
                    % Initialize the carrier resource grid for one slot and map PRS symbols onto
                    % the grid
                    txSlotGrid = nrResourceGrid(cell_carrier{gNBIdx});
                    txSlotGrid(prsInd{gNBIdx}{resIdx}) = db2mag(powerPRS)*prsSym{gNBIdx}{resIdx};
                    reshapedSymb = reshape(txSlotGrid,[],ports);
                    
                    % Get the transmit beam index
                    beamIdx = find(activeRes == resIdx);
                    
                    % Apply the digital beamforming
                    if ~isempty(beamIdx)
                        bfSymb = reshapedSymb * wT{gNBIdx}(:,beamIdx)';
                        bfGrid = bfGrid + reshape(bfSymb,size(bfGrid));
                    end
                end
                
                % Perform OFDM Modulation
                % Generate the time-domain waveform by performing the OFDM modulation.
                
                tbfWaveform = nrOFDMModulate(cell_carrier{gNBIdx},bfGrid);
                
                % Normalize the beamformed time-domain waveform over the number of transmit
                % antennas
                
                tbfWaveform = tbfWaveform/sqrt(numTxAnt);
                
                % Send the Waveform through the Channel
                % Append zeros at the end of the transmitted waveform to flush the channel 
                % content and then pass the time-domain waveform through the scattering MIMO channel. 
                % These zeros take into account any delay introduced in the channel.
                
                % Append zeros to the transmit waveform to account for channel delay
                tbfWaveform = [tbfWaveform; zeros(maxChDelay,numTxAnt)];
                % Pass the waveform through the channel
                fadWave = cell_channel{gNBIdx}(tbfWaveform);
                
                fadWaveG = fadWave*min(rxGain);
                noise = N0DL*complex(randn(size(fadWaveG)),randn(size(fadWaveG)));
                rxWaveform = fadWaveG + noise;
%                     rxWaveform = fadWave;   
                
                % Timing Synchronization
                % Perform the timing synchronization by cross correlating the 
                % received reference symbols with a local copy of PRS symbols.
                
                % Generate reference symbols and indices
                refSym = nrPRS(cell_carrier{gNBIdx},prs(gNBIdx));
                refInd = nrPRSIndices(cell_carrier{gNBIdx},prs(gNBIdx));
                
                % Estimate timing offset
                offset = nrTimingEstimate(cell_carrier{gNBIdx},rxWaveform,refInd,refSym);
                if offset > maxChDelay
                    offset = 0;
                end
                % Correct timing offset
                syncTdWaveform = rxWaveform(1+offset:end,:);
                
                % OFDM Demodulation and Receive Beamforming
                % OFDM Demodulation
                % OFDM demodulate the synchronized time-domain waveform.
                
                rxGrid = nrOFDMDemodulate(cell_carrier{gNBIdx},syncTdWaveform);
                
                % Calculate the steering vector for the angle of reception.
                
                wR{gNBIdx} = SteerVecRx(fc,coarse_AoA(gNBIdx,:)');
                
                temp = rxGrid;
                if strcmpi(prm.FreqRange,'FR1')
                    % Beamforming without combining
                    rbfGrid = reshape(reshape(temp,[],numRxAnt).*wR{gNBIdx}',size(temp,1),size(temp,2),[]);
                else % 'FR2'
                    % Beamforming with combining
                    rbfGrid = reshape(reshape(temp,[],numRxAnt)*conj(wR{gNBIdx}),size(temp,1),size(temp,2),[]);
                end
                %% Refined Beam Determination
                % After the OFDM demodulation, the UE measures the RSRP for all the PRS 
                % resources transmitted in different beams, given the current receive beam. 
                % Perform these measurements by using the helper function hPRSMeasurements.
                
                % Perform RSRP measurements
                meas = hPRSMeasurements(cell_carrier{gNBIdx},prs(gNBIdx),rbfGrid);
                
                % Display the measurement quantities for all PRS resources in dBm
                RSRPdBm = meas.RSRPdBm;
                if plot_rays
                    disp(['RSRP measurements (in dBm) of gNB ' num2str(reachable_gNBs(gNBIdx).pci) ':' 13 num2str(RSRPdBm)]);
                end
                % Identify the maximum RSRP value from the measurements and find the best corresponding beam.
                
                % Get the transmit beam index with maximum RSRP value
                [~,maxRSRPIdx] = max(RSRPdBm(logical(prsTransmitted)));
                
                % Get the PRS resource index with maximum RSRP value
                [~,maxRSRPResIdx] = max(RSRPdBm);
                if numBeams == 0
                    disp('Refinement has not happened, as PRS is not transmitted')
                else
                    refBeamWts = wT{gNBIdx}(:,maxRSRPIdx);
                    [prsAzBeamWidth, AzRange{gNBIdx}] = beamwidth(reachable_gNBs(gNBIdx).antennaArray,fc,'PropagationSpeed',c,'Weights',refBeamWts,'Cut','Elevation','CutAngle',prsBeamAng(1,maxRSRPIdx));
                    prsrsElBeamWidth = beamwidth(reachable_gNBs(gNBIdx).antennaArray,fc,'PropagationSpeed',c,'Weights',refBeamWts,'Cut','Azimuth','CutAngle',prsBeamAng(2,maxRSRPIdx));
                    
                    refined_AoD(gNBIdx, :) = prsBeamAng(:,maxRSRPIdx);
                    refined_AoA(gNBIdx, :) = refined_AoD(gNBIdx, :)+reachable_gNBs(gNBIdx).orientation;
                    refined_AoA(gNBIdx, 1) = refined_AoA(gNBIdx, 1)+180;
                    if (refined_AoA(gNBIdx, 1)>180)
                        refined_AoA(gNBIdx, 1) = refined_AoA(gNBIdx, 1)-360;
                    end
                    
                    refined_AoA(gNBIdx,2) = -refined_AoA(gNBIdx, 2);
                    rho_aoa = refined_AoA.*pi/180;
                    [~, real_angle] = rangeangle(reachable_gNBs(gNBIdx).getXYZ()', UE_xyz');
                    if plot_rays
                        disp([ 'gNB ' num2str(reachable_gNBs(gNBIdx).pci) ' with transmit-end beam refinement:' 13 'Refined transmit beam ('...
                            num2str(maxRSRPIdx) ') corresponds to PRS resource '...
                            num2str(maxRSRPResIdx) ' is selected in the direction ['...
                            num2str(refined_AoD(gNBIdx,1)) ';' num2str(refined_AoD(gNBIdx,2))...
                            ']' 13 'after the selected SSB in direction [' ...
                            num2str(coarse_AoD(gNBIdx,1)) ';' num2str(coarse_AoD(gNBIdx,2))...
                            ']' 13 'the refined AoA direction is [' ...
                            num2str(refined_AoA(gNBIdx,1)) ';' num2str(refined_AoA(gNBIdx,2))...
                            ']' 13 'the real AoA direction is [' ...
                            num2str(real_angle(1)) ';' num2str(real_angle(2))...
                            ']']);
                    end
                end
            else
                coarse_AoA(gNBIdx,:) =coarse_AoD(gNBIdx, :)+reachable_gNBs(gNBIdx).orientation;
                coarse_AoA(gNBIdx,1) = coarse_AoA(gNBIdx, 1)+180;
                if (coarse_AoA(gNBIdx, 1)>180)
                    coarse_AoA(gNBIdx, 1) = coarse_AoA(gNBIdx, 1)-360;
                end
                coarse_AoA(gNBIdx,2) = -coarse_AoA(gNBIdx, 2);
                rho_aoa = coarse_AoA.*pi/180;
            end
        end
        
        estimated_ue = Non_linear_LS_AoA2(rho_aoa, u_0, s_aoa, K, stop_cond, dimensions)
    end
end
