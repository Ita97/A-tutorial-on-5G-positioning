clc;  close all;
clearvars -except viewer;

%% Simulation Variables
rng(42);                       % Set RNG state for repeatability
pd = truncate(makedist("Normal", "sigma",50),-100,100); % probability distribution
                            
reflections_order = 2;       % 0 (LOS) default | 1 (NLOS)
diffractions_order = 0;      % 0 default | 1 | 2
angular_sep = "medium";         % "medium" default | "high" | "low"
plot_rays = true;            % plot raytracer rays on map
refinement = false;          % true (SPI x ToA and P2 x AoA) default | false

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
TxArraySize = [4 4];        % BS Antenna Array Size - [vertical, horizontal]
RxArraySize = [2 2];        % UE Antenna Array Size - [vertical, horizontal]
numTxAnt = prod(TxArraySize);
numRxAnt = prod(RxArraySize);
TxAzimuthRange = [-60 60];
TxElevationRange = [-90 0];
RxAzimuthRange = [-180 180];
RxElevationRange = [0 90];
%ueArrayOrientation = [0 45].';      % azimuth (0 deg is E, 90 deg is N) and 
                                    % elevation (positive points upwards) in deg

ueArrayOrientation = [-60 0].';

%% UE positions
load ue_positions/giuriati.mat
%listUEPos = [45.479654,9.230301];
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
    ueSite{UEIdx}.Antenna = ueArray{UEIdx};
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

%% For each UE position
for posIdx = 1: numUEPos
    UEPos = listUEPos(posIdx,:);
    UE_xyz = list_UE_xyz(posIdx, :);

    show(ueSite{posIdx});

    %% Ray Tracing Analysis
    % Perform ray tracing analysis using the Shooting and Bouncing Rays (SBR) method. 
    % The SBR method includes effects from surface reflections and diffractions but 
    % does not include effects from refraction or scattering.
    pm = propagationModel("raytracing", "Method", "sbr","MaxNumReflections",reflections_order,"MaxNumDiffractions",diffractions_order,"AngularSeparation",angular_sep);
    [channel, reachable_gNBs, reachable_gNBs_idx] = CDLChannelConfig(gNBs, bsSite, ueSite{posIdx},ueSite{ref_site}, pm, ofdmInfo, plot_rays);
    numgNBs = numel(reachable_gNBs);
      
    los = false(numgNBs, 1);
    for i= 1: numgNBs
        los(i) = channel{reachable_gNBs_idx(i,1), reachable_gNBs_idx(i,2)}.HasLOSCluster;
    end
    %% Drop NLOS gNBs
    if drop_nlos
        if any(not(los))
            reachable_gNBs_idx = reachable_gNBs_idx(los,:);
            reachable_gNBs = reachable_gNBs(los);
            %numgNBs = numel(reachable_gNBs);
        end
    end

    if numel(reachable_gNBs)< max_gNBs 
        numgNBs = numel(reachable_gNBs);  
    else 
        numgNBs=max_gNBs; 
        reachable_gNBs_idx = reachable_gNBs_idx(1:numgNBs,:);
        reachable_gNBs = reachable_gNBs(1:numgNBs);
    end  

    disp([13 'For position ' num2str(posIdx) ' the number of available gNB cells are: ' num2str(numgNBs)])

    %% PRS Configuration
    % Configure PRS parameters for all of the gNBs. 
    % Configure the PRS for different gNBs such that no overlap exists to avoid the problem of hearability. 
    % This example considers different slot offsets for the PRS from different gNBs to avoid the overlap among PRS signals. 
    % This overlapping can be avoided in multiple ways by configuring the time-frequency aspects of PRS signals in an appropriate way 
    % (like choosing the nonoverlapped symbol allocation or frequency allocation, or by using muting pattern configuration parameters). 
    
    
    % Slot offsets of different PRS signals
    prsSlotOffsets = 0:2:(2*numgNBs - 1);
    prsIDs = randperm(4096,numgNBs) - 1;
    
    % Configure PRS properties
    prs = nrPRSConfig;
    prs.PRSResourceSetPeriod = [numgNBs*2 0];
    prs.PRSResourceOffset = 0;
    prs.PRSResourceRepetition = 1;
    prs.PRSResourceTimeGap = 1;
    prs.MutingPattern1 = [];
    prs.MutingPattern2 = [];
    prs.NumRB = NSizeGrid;
    prs.RBOffset = 0;
    prs.CombSize = 12;
    prs.NumPRSSymbols = 12;
    prs.SymbolStart = 0;
    prs = repmat(prs,1,numgNBs);
    for gNBIdx = 1:numgNBs
        prs(gNBIdx).PRSResourceOffset = prsSlotOffsets(gNBIdx);
        prs(gNBIdx).NPRSID = prsIDs(gNBIdx);
    end
    
    %% PDSCH Configuration

    % It assumes that the data is transmitted from all of the gNBs that have a single transmission layer.
    pdsch = nrPDSCHConfig;
    pdsch.PRBSet = 0:51;
    pdsch.SymbolAllocation = [0 carrier(1,1).SymbolsPerSlot];
    pdsch.DMRS.NumCDMGroupsWithoutData = 1;
    pdsch = repmat(pdsch,1,numgNBs);
    
    % validate the number of layers
    numLayers = [pdsch(:).NumLayers];
    if ~all(numLayers == 1)
        error('nr5g:invalidNLayers',['The number of transmission layers ' ...
            'configured for the data transmission must be 1.']);
    end

    % Resource Generation
    cell_carrier = cell(numgNBs,1);
    for i = 1:numgNBs
            cell_carrier{i} = carrier(reachable_gNBs_idx(i,1), reachable_gNBs_idx(i,2));
    end
    
    [prsGrid, DLdataGrid] = PRSResourceGeneration(prs, pdsch, cell_carrier, numTxAnt, nFrames);

%% For each Test
    for test = 1:num_tests

        disp(['Test #' num2str(test)]);
        if length(reachable_gNBs) < max_gNBs 
            numgNBs = length(reachable_gNBs);  
        else 
            numgNBs=max_gNBs; 
        end
        u_0 = randi([-100,100],1,dimensions);
        if dimensions==3
            u_0(1,3) = 0;
        end
        
        %% Add Signal Delays
    
        sampleDelay = cell(1,numgNBs);
        cell_channel = cell(1,numgNBs);
        for i = 1:numgNBs
            if not(BS_synch)
                desynch = random(pd)*1e-9;
            else
                desynch = 0;
            end
            cell_channel{i} = channel{reachable_gNBs_idx(i,1),reachable_gNBs_idx(i,2)};
            chInfo = info(cell_channel{i});
            sampleDelay{i} = ceil((chInfo.PathDelays+desynch)*cell_channel{i}.SampleRate);
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

        %% Perform PRS OFDM Modulation
        % Perform OFDM modulation of PRS and data signal at each gNB
        txWaveform = cell(1,numgNBs);
        for gNBIdx = 1:numgNBs
            release(cell_channel{gNBIdx});
            if cell_channel{gNBIdx}.TransmitAndReceiveSwapped
                swapTransmitAndReceive(cell_channel{gNBIdx});
            end
            cell_channel{gNBIdx}.RandomStream = 'Global stream';
            cell_carrier{gNBIdx}.NSlot = 0;
            txWaveform{gNBIdx} = nrOFDMModulate(cell_carrier{gNBIdx},prsGrid{gNBIdx} + DLdataGrid{gNBIdx});
        end

        %% Signal through Channel - Downlink
        % Model the received signal at the UE by delaying each gNB 
        % transmission according to the values in sampleDelay and by attenuating 
        % the received signal from each gNB. Ensure that all of the waveforms 
        % are of the same length by padding the received waveform from each gNB 
        % with relevant number of zeros.

        rxWaveform = signal2Channel(true, TxPower, txWaveform, cell_channel, sampleDelay, numgNBs, los, N0DL, ofdmInfo.Nfft, NSizeGrid);

        %% DL-TOA
        % The UE performs the correlation on the incoming signal with reference PRS generated for each gNB, 
        % and cellsToBeDetected number of best cells are selected based on correlation outcome. 
        % Compute TOAs of the signals from each gNB by using nrTimingEstimate function.         
        [rho_dltoa, s_dltoa] = toa_estimation(true, rxWaveform, prsGrid, reachable_gNBs, reachable_gNBs_idx, cell_carrier, ofdmInfo, reflections_order, refinement, true);
        
        %% DL-TDOA
        % To estimate the position of UE, compute RSTD values by using the absolute TOAs corresponding to the best cells.
        [rho_tdoa, s_tdoa, ~, ~]=compute_tdoa_measurements2(s_dltoa, rho_dltoa);
        estimated_ue = Non_linear_LS_TDOA(rho_tdoa, u_0, s_tdoa, K, stop_cond, dimensions)
        
        %% Bidimensional Error
        error = sqrt(sum((estimated_ue(1:2)-UE_xyz(1:2)).^2))

    end

end


