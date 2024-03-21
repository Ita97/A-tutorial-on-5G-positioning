function [channel, reachable_gNBs, reachable_gNBs_idx] = CDLChannelConfig(gNBs, bsSite, ueSite, refUeSite, pm, ofdmInfo, plot_rays)
    numgNBs = numel(gNBs);
    n_cell = numel(gNBs{1});
    channel = cell(numgNBs,n_cell);
    scatPos=cell(numgNBs,n_cell);
    LOSRay = cell(numgNBs,n_cell);
    los = zeros(1, numgNBs);
    reachable_gNBs_idx = [];
    reachable_gNBs = [];
    
    for gNBIdx = 1:numgNBs
        rays = raytrace(bsSite{gNBIdx},ueSite,pm,"Type","pathloss");
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
                scatterPos_xyz = from_lat_long_to_xyz(SPSite{i}, refUeSite);
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
            channel{gNBIdx, cellIdx}.CarrierFrequency = bsSite{gNBIdx}.TransmitterFrequency;
            channel{gNBIdx, cellIdx}.NormalizeChannelOutputs = false; % do not normalize by the number of receive antennas, this would change the receive power
            channel{gNBIdx, cellIdx}.NormalizePathGains = false;      % set to false to retain the path gains
            
            channel{gNBIdx, cellIdx}.ReceiveAntennaArray = ueSite.Antenna;
            channel{gNBIdx, cellIdx}.ReceiveArrayOrientation = [ueSite.AntennaAngle(1); (-1)*ueSite.AntennaAngle(2); 0];  % the (-1) converts elevation to downtilt
            
            channel{gNBIdx, cellIdx}.TransmitAntennaArray = bsSite{gNBIdx}(cellIdx).Antenna;
            channel{gNBIdx, cellIdx}.TransmitArrayOrientation = [bsSite{gNBIdx}(cellIdx).AntennaAngle(1); (-1)*bsSite{gNBIdx}(cellIdx).AntennaAngle(2); 0];   % the (-1) converts elevation to downtilt 
        
            channel{gNBIdx, cellIdx}.SampleRate = ofdmInfo.SampleRate;
        end
    
    
    end
        
    if numel(logical(los)) ~= numgNBs
            error('nr5g:InvalidLOSLength',['Length of line of sight flag (' num2str(numel(logical(los))) ...
            ') must be equal to the number of configured gNBs (' num2str(numgNBs) ').']);
    end

end