function meas = hPRSMeasurements(carrier,prs,grid)

% 
%hPRSMeasurements PRS based reference signal measurements
%   MEAS = hCSIRSMeasurements(CARRIER,PRS,GRID) returns the PRS based
%   reference signal measurements (reference signal received power (RSRP),
%   received signal strength indicator (RSSI) and reference signal received
%   quality (RSRQ)) in a structure MEAS, as described in
%
%   CARRIER is a carrier specific configuration object as described in
%   <a href="matlab:help('nrCarrierConfig')">nrCarrierConfig</a> with the following properties:
%
%   NCellID           - Physical layer cell identity
%   SubcarrierSpacing - Subcarrier spacing in kHz
%   CyclicPrefix      - Cyclic prefix
%   NSizeGrid         - Number of resource blocks (RBs) in carrier resource
%                       grid
%   NStartGrid        - Start of carrier resource grid relative to common
%                       resource block 0 (CRB 0)
%   NSlot             - Slot number
%   NFrame            - System frame number
%
%   PRS is a PRS specific configuration object as described in
%   <a href="matlab:help('nrPRSConfig')">nrPRSConfig</a> with the following properties specified for one or more
%   PRS resources:
%
%   NumRB                 - Number of RBs allocated for a PRS resource
%   RBOffset              - Starting RB index of CSI-RS allocation
%                           relative to carrier resource grid
%   NID                   - Scrambling identity
%
%   GRID is a 3-dimensional array of the resource elements, for one slot
%   across all receive antennas. GRID is of size K-by-L-by-R array, where K
%   represents the number of subcarriers, L represents the number of OFDM
%   symbols, and R represents the number of receive antennas.
%
%   MEAS is a structure including the fields:
%   RSRPPerAntennaPerResource - Matrix of linear RSRP values with rows
%                               corresponding to receive antennas and
%                               columns corresponding to CSI-RS resources
%   RSSIPerAntennaPerResource - Matrix of linear RSSI values with rows
%                               corresponding to receive antennas and
%                               columns corresponding to CSI-RS resources
%   RSRQPerAntennaPerResource - Matrix of linear RSRQ values with rows
%                               corresponding to receive antennas and
%                               columns corresponding to CSI-RS resources
%   RSRP                      - Row vector of reported RSRP values for all
%                               CSI-RS resources, where each element
%                               represents the maximum value of linear
%                               RSRPs of all receive antennas
%   RSSI                      - Row vector of reported RSSI values for all
%                               CSI-RS resources, where each element
%                               represents the RSSI value of the receive
%                               antenna with the maximum RSRQ
%   RSRQ                      - Row vector of reported RSRQ values for all
%                               CSI-RS resources, where each element
%                               represents the maximum value of linear
%                               RSRQs of all receive antennas
%   RSRPdBm                   - The RSRP expressed in deciBels relative to
%                               1 milliwatt (in dBm)
%   RSSIdBm                   - The RSSI expressed in deciBels relative to
%                               1 milliwatt (in dBm)
%   RSRQdB                    - The RSRQ expressed in deciBels (in dB)
%
%   Note that the columns of above fields are in the same order, as the
%   resources are specified in the input configuration object CSIRS.
%

    % Get the number of PRS resources

    numPRSRes = max([numel(prs.NumPRSSymbols), numel(prs.PRSResourceOffset), ...
        numel(prs.SymbolStart), numel(prs.REOffset), numel(prs.NPRSID)]);
    SymbolLocations = 1;

    % Number of receive antennas
    nRx = size(grid,3);

    % Initialize the variables for measurements
    rsrpPerAntPerRes = zeros(nRx,numPRSRes); % Row indicates the receive antenna index and column indicates the CSI-RS resource index
%     rssiPerAntPerRes = zeros(nRx,numPRSRes);
%     rsrqPerAntPerRes = zeros(nRx,numPRSRes);
    
    % Reference indices and symbols generation
    
    tmpRefSym = nrPRS(carrier,prs, "OutputResourceFormat","cell");
    refIndLin = nrPRSIndices(carrier,prs, "OutputResourceFormat","cell");
    
    for rxAntIdx = 1:nRx % Loop over all receive antennas
        gridRxAnt = grid(:,:,rxAntIdx);
        for resIdx = 1:numPRSRes % Loop over all PRS resources for one receive antenna
            if isscalar(prs.NumRB)
                N = prs.NumRB;
            else
                N = prs.NumRB(resIdx);
            end

            if isscalar(prs.RBOffset)
                rbOffset = prs.RBOffset;
            else
                rbOffset = prs.RBOffset(resIdx);
            end
            if iscell(SymbolLocations)
                l0 = SymbolLocations{resIdx}(1); % l0 value
            else
                l0 = SymbolLocations(1);
            end
            % For RSSI measurement, generate the indices of all resource elements in OFDM symbols containing PRS resource 
%             prsSymIndices = l0 + info.LPrime{resIdx} + 1; % 1-based symbol indices
%             numPRSSym = numel(prsSymIndices);
%             rssiIndices = repmat((1:N*12).' + rbOffset*12,1,numPRSSym) + repmat((prsSymIndices - 1)*12*carrier.NSizeGrid,12*N,1); % 1-based, linear carrier-oriented indices
%             % Extract the modulation symbols using the indices, which
%             % corresponds to RSSI measurement
%             rssiSym = gridRxAnt(rssiIndices);

            numPorts = 1;
            if numPorts > 2
                ports = 2; 
                tmp2 = reshape(tmpRefSym{resIdx},[],numPorts);
                refSym = reshape(tmp2(:,1:2),[],1);
            else
                ports = numPorts; 
                refSym = tmpRefSym{resIdx};
            end
            indTmp = reshape(refIndLin{resIdx},[],numPorts);
            indTmp = double(indTmp(:,1:ports)) - (0:ports-1)*carrier.NSizeGrid*12*carrier.SymbolsPerSlot; % Port 3000 or Ports 3000 and 3001
            % Extract the received CSI-RS symbols using locally generated
            % CSI-RS indices
            rxSym = reshape(gridRxAnt(indTmp),[],1);

            if ~isempty(refSym)
                % Perform the following calculations for NZP-CSI-RS
                % resources only
                rsrpPerAntPerRes(rxAntIdx,resIdx) = abs(mean(rxSym.*conj(refSym))*ports)^2;
%                 rssiPerAntPerRes(rxAntIdx,resIdx) = sum(rssiSym.*conj(rssiSym))/numPRSSym;
%                 rsrqPerAntPerRes(rxAntIdx,resIdx) = N*rsrpPerAntPerRes(rxAntIdx,resIdx)/rssiPerAntPerRes(rxAntIdx,resIdx);
            end
        end
    end

    rsrpPerAntPerRes(isnan(rsrpPerAntPerRes)) = 0;
%     rssiPerAntPerRes(isnan(rssiPerAntPerRes)) = 0;
%     rsrqPerAntPerRes(isnan(rsrqPerAntPerRes)) = 0;

    % Extract the maximum value of per antenna measurements (values
    % correspond to individual receiver branches) for all CSI-RS resources
    % to get the RSRP, RSRQ reporting values, as described in TS 38.215
    % Sections 5.1.2 and 5.1.4
    rsrp = max(rsrpPerAntPerRes,[],1);
%     [rsrq,rsrqIdx] = max(rsrqPerAntPerRes,[],1);

    % Extract RSSI value of the receive antenna with the maximum RSRQ for
    % all CSI-RS resources
%     rssi = rssiPerAntPerRes(rsrqIdx + (0:numPRSRes-1)*nRx);

    % Convert the RSRP, RSSI and RSRQ reporting values from linear scale to
    % logarithmic scale
    RSRPdBm = 10*log10(rsrp) + 30; % dBm
%     RSSIdBm = 10*log10(rssi) + 30; % dBm
%     RSRQdB = 10*log10(rsrq);       % dB

    % Create an output structure with the fields representing the
    % measurements
    meas.RSRPPerAntennaPerResource = rsrpPerAntPerRes;
%     meas.RSSIPerAntennaPerResource = rssiPerAntPerRes;
%     meas.RSRQPerAntennaPerResource = rsrqPerAntPerRes;
    meas.RSRP = rsrp;
%     meas.RSSI = rssi;
%     meas.RSRQ = rsrq;
    meas.RSRPdBm = RSRPdBm;
%     meas.RSSIdBm = RSSIdBm;
%     meas.RSRQdB = RSRQdB;



end