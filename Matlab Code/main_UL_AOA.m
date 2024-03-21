clc;  close all;
clearvars -except viewer;

%% Simulation Variables

rng(42);                       % Set RNG state for repeatability
pd = truncate(makedist("Normal", "sigma",50),-100,100); % probability distribution
                            
reflections_order = 2;       % 0 (LOS) default | 1 (NLOS)
plot_rays = true;            % plot raytracer rays on map
refinement = false;          % true (parabolic interpolation) default | false

num_tests = 1;             
dimensions = 3;            % 2/3-D
K = 1000;                  % number of iteration for NLS algorithm
stop_cond = 1e-4;

drop_nlos = true;           % drop NLOS BS
BS_synch = true;            % BS synchronized?
max_gNBs = 4;               % Max number of BS used

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
BSTxPower = 33;               % BS Tx Power [dBm] Urban
%BSTxPower = 49;               % BS Tx Power [dBm] Rural
UETxPower = 23;               % UE Tx Power [dBm]

ref_site = 1;               % BS used as reference for geo<->cartesian conversion
TxArraySize = [8 8];        % BS Antenna Array Size
RxArraySize = [2 2];        % UE Antenna Array Size 
numBSAnt = prod(TxArraySize);
numUEAnt = prod(RxArraySize);
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

for gNBIdx = 1 : numgNBs 
    bsSite{gNBIdx} = txsite("Name","BS"+gNBIdx, ...
        "Latitude",gNBs{gNBIdx}(:).getCoordinates(1),"Longitude",gNBs{gNBIdx}(:).getCoordinates(2),...
        "AntennaAngle",gNBs{gNBIdx}(:).getOrientation()',...
        "AntennaHeight",h_gNB,...  % in m
        "TransmitterFrequency",fc, ...
        "TransmitterPower",db2pow(BSTxPower-30),...
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
ue_carrier = nrCarrierConfig;
% Configure carrier properties
ue_carrier.NSizeGrid = NSizeGrid;
ue_carrier.SubcarrierSpacing = SCS;
validateCarriers(ue_carrier);
ofdmInfo = nrOFDMInfo(ue_carrier);

%% Noise

NF = 10^(NoiseFigureRx/10);
Te = T_ant + 290*(NF-1);    % noise temperature [K]
kB = physconst("Boltzmann"); % [J/K]
N0DL = sqrt((kB*ofdmInfo.SampleRate*Te)/2);     % noise amplitude DL
%N0DL=0;
NF = 10^(NoiseFigureTx/10);
Te = T_ant + 290*(NF-1);    % noise temperature [K]
N0UL = sqrt((kB*ofdmInfo.SampleRate*Te)/2);     % noise amplitude UL

% Shadowing Fading (SF)
sigma_SF = 3; % dB

%% For each UE position
for posIdx = 1: numUEPos
    numgNBs = size(gNBs,1);
    UEPos = listUEPos(posIdx,:);
    UE_xyz = list_UE_xyz(posIdx, :); 

    show(ueSite{posIdx});

    %% Ray Tracing Analysis
    % Perform ray tracing analysis using the Shooting and Bouncing Rays (SBR) method. 
    % The SBR method includes effects from surface reflections and diffractions but 
    % does not include effects from refraction or scattering.
    pm = propagationModel("raytracing", "Method", "sbr","MaxNumReflections",reflections_order);
    
    % Channel Modelling
    [channel, reachable_gNBs, reachable_gNBs_idx] = CDLChannelConfig(gNBs, bsSite, ueSite{posIdx},ueSite{ref_site}, pm, ofdmInfo, plot_rays);
    numgNBs = numel(reachable_gNBs);        % we consider from now on, only the reachable cell

    los = false(numgNBs, 1);
    for i= 1: numgNBs
        los(i) = channel{reachable_gNBs_idx(i,1), reachable_gNBs_idx(i,2)}.HasLOSCluster;
    end
    %% Drop NLOS gNBs
    if drop_nlos
        % LOS = false(numgNBs, 1);
        % for i= 1: numgNBs
        %     LOS(i) = channel{reachable_gNBs_idx(i,1), reachable_gNBs_idx(i,2)}.HasLOSCluster;
        % end
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

    %% SRS configuration
    
    srs = nrSRSConfig;
    srs.NumSRSSymbols = 8;          % Number of OFDM symbols allocated per slot (1, 2, 4, 8, 12)
    srs.SymbolStart = 0;            % Starting OFDM symbol within a slot
    srs.NumSRSPorts = 1;            % Number of SRS antenna ports (1,2,4).
    srs.SRSPositioning = 1;         % Rel-16
    srs.FrequencyStart = 0;         % Frequency position of the SRS in BWP in RBs
    srs.NRRC = 0;                   % Additional offset from FreqStart specified in blocks of 4 PRBs (0...67)
    srs.CSRS = 63;                  % Bandwidth configuration C_SRS (0...63). It controls the allocated bandwidth to the SRS.. CHECK THE TABLE!
    srs.BSRS = 0;                   % Bandwidth configuration B_SRS (0...3). It controls the allocated bandwidth to the SRS
    srs.BHop = 0;                   % Frequency hopping configuration (0...3). Set BHop < BSRS to enable frequency hopping
    srs.KTC = 8;                    % Comb number (2,4,8). Frequency density in subcarriers
    srs.Repetition = 2;             % Repetition (1,2,4). It disables frequency hopping in blocks of |Repetition| symbols
    srs.SRSPeriod = [2 0];          % Periodicity and offset in slots. SRSPeriod(2) must be < SRSPeriod(1)
    srs.ResourceType = 'periodic';  % Resource type ('periodic', 'semi-persistent','aperiodic'). Use 'aperiodic' to disable inter-slot frequency hopping
    
    %% PUSCH Configuration
    % It assumes that the data is transmitted from all of the gNBs that have a single transmission layer.
    pusch = nrPUSCHConfig;
    pusch.PRBSet = 0:51;
    pusch.SymbolAllocation = [0 ue_carrier.SymbolsPerSlot];
    pusch.DMRS.NumCDMGroupsWithoutData = 1;
    
    % validate the number of layers
    numLayers = [pusch(:).NumLayers];
    if ~all(numLayers == 1)
        error('nr5g:invalidNLayers',['The number of transmission layers ' ...
            'configured for the data transmission must be 1.']);
    end

    % Precoding
    cell_channel = cell(1,numgNBs);
    newWtx = cell(1,numgNBs);
    for gNBIdx = 1:numgNBs 
        cell_channel{gNBIdx} = channel{reachable_gNBs_idx(gNBIdx,1),reachable_gNBs_idx(gNBIdx,2)};
        % configure for UL transmission assuming channel reciprocity
        swapTransmitAndReceive(cell_channel{gNBIdx});
        estChannelGrid = getInitialChannelEstimate(ue_carrier,numUEAnt,cell_channel{gNBIdx});
        newWtx{gNBIdx} = getPrecodingMatrix(ue_carrier,pusch,estChannelGrid,[]); 
        % PUSCH PRB bundling (TS 38.214 Section 6.1.2.3), 2, 4, or [] to signify "wideband"   
    end

    % Resource Generation
    [srsGrid, ULdataGrid] = SRSResourceGeneration(srs, pusch,ue_carrier, numUEAnt, nFrames, newWtx);
    
    %% For each Test
    for test = 1:num_tests

        disp(['Test #' num2str(test)]);
        numgNBs = numel(reachable_gNBs);
        u_0 = randi([-100,100],1,dimensions);
        if dimensions==3
            u_0(1,3) = 0;
        end
        
        %% Add Signal Delays
        sampleDelay = zeros(1,numgNBs);
        for i = 1:numgNBs
            chInfo = info(cell_channel{i});
            sampleDelay(i) = ceil(max(chInfo.PathDelays*cell_channel{i}.SampleRate));
        end
        maxChDelay = max(sampleDelay);

        if test>1 % return back the channel order
            cell_carrier = cell(numgNBs,1);
            cell_channel = cell(numgNBs,1);
            for i = 1:numgNBs
                idx = [ reachable_gNBs_idx(i,1), reachable_gNBs_idx(i,2) ];
                cell_carrier{i} = carrier(idx(1), idx(2));
                cell_channel{i} = channel{idx(1), idx(2)};
                %reset(cell_channel{i});
            end
        end

        %% Uplink Measurements Extraction
        
        %% Perform SRS OFDM Modulation
        % Perform OFDM modulation of SRS and data signal at each gNB
        txWaveform = cell(1,numgNBs);
        for gNBIdx = 1:numgNBs
            release(cell_channel{gNBIdx});
            cell_channel{gNBIdx}.RandomStream = 'Global stream';
            ue_carrier.NSlot = 0;
            txWaveform{gNBIdx} = nrOFDMModulate(ue_carrier,srsGrid{gNBIdx} + ULdataGrid{gNBIdx});
        end

        %% Signal through Channel - Uplink
        % Model the received signal at each gNB, delaying the 
        % transmission according to the values in sampleDelay and by attenuating 
        % the received signal of each gNB. 
        %TxAmp = db2mag(UETxPower-30)*sqrt(ofdmInfo.Nfft^2/(NSizeGrid*12*numUEAnt));
        %rxWaveform = signal2Channel(false, TxAmp, txWaveform, cell_channel, sampleDelay, numgNBs, numBSAnt, numUEAnt, los, N0UL);
        rxWaveform = cell(1,numgNBs);
        for gNBIdx = 1: numgNBs
            txWaveform{gNBIdx} = [txWaveform{gNBIdx}; zeros(sampleDelay(gNBIdx),size(txWaveform{gNBIdx},2))];
            [rxWaveform{gNBIdx},~,~] = cell_channel{gNBIdx}(txWaveform{gNBIdx});
            spLoss = -max(cell_channel{gNBIdx}.AveragePathGains);
    
            % add noise
            rxGain = 10.^(spLoss/20);
            noise = N0UL*complex(randn(size(rxWaveform{gNBIdx})),randn(size(rxWaveform{gNBIdx})));
            rxWaveform{gNBIdx} = rxWaveform{gNBIdx}*min(rxGain)+noise;
        
            % Reset to DL configuration
            swapTransmitAndReceive(cell_channel{gNBIdx});
        end
        
        %% UL-AOA
        [rho_aoa, s_aoa] = compute_aoa_measurements(rxWaveform, reachable_gNBs)
        estimated_ue = Non_linear_LS_AoA2(rho_aoa, u_0, s_aoa, K, stop_cond, dimensions)

    
    end

end