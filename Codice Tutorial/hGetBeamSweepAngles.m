function beamAng = hGetBeamSweepAngles(numBeams,azSweepRange,elSweepRange,azBW,elBW,varargin)
%hGetBeamSweepAngles Return azimuth and elevation angle pairs for all beams
%
%   BEAMANG = hGetBeamSweepAngles(NUMBEAMS,AZSWEEPRANGE,ELSWEEPRANGE, ...
%   AZBW,ELBW) returns azimuth and elevation angle pairs BEAMANG, given the
%   number of beams NUMBEAMS, angular sweep range in azimuth plane
%   AZSWEEPRANGE, angular sweep range in elevation plane ELSWEEPRANGE,
%   beamwidth of antenna array in azimuth plane AZBW, and beamwidth of
%   antenna array in elevation plane ELBW.
%   The beams are equally spaced within the sweeps in both azimuth and
%   elevation directions.
%
%   BEAMANG = hGetBeamSweepAngles(...,ELSWEEP) also allows a boolean flag
%   ELSWEEP to enable or disable the elevation sweep.

%   Copyright 2020 The MathWorks, Inc.

    narginchk(5,6);
    if nargin == 5
        elSweep = true;
    else
        elSweep = varargin{1};
    end
    if ~elSweep
        elSweepRange = [0 0];
    end

    if numBeams == 0
        beamAng = [];
    else
        if numel(unique(azSweepRange)) == 2 && numel(unique(elSweepRange)) == 2
            % Both azimuth and elevation sweep
            effNumAzBeams = ceil((azSweepRange(2) - azSweepRange(1))/azBW);
            effNumElBeams = ceil((elSweepRange(2) - elSweepRange(1))/elBW);
            effAzElBeamsRatio = effNumAzBeams/effNumElBeams;

            temp = numBeams./(numBeams:-1:1);
            divVals = temp(temp == floor(temp))';
            pairsAzEl = [divVals flipud(divVals)];
            [~,index] = min(abs((pairsAzEl(:,1)./pairsAzEl(:,2)) - effAzElBeamsRatio));

            beamsInAz = pairsAzEl(index,1);
            if beamsInAz == 1
                azDir = (azSweepRange(1) + azSweepRange(2))/2;
            else
                azBeamSeparation = (azSweepRange(2) - azSweepRange(1))/beamsInAz;
                azDir = azSweepRange(1) + azBeamSeparation:azBeamSeparation:azSweepRange(2);
            end
            
            beamsInEl = pairsAzEl(index,2);
            if beamsInEl == 1
                elDir = (elSweepRange(1) + elSweepRange(2))/2;
            else
                elBeamSeparation = (elSweepRange(2) - elSweepRange(1))/beamsInEl;
                elDir = elSweepRange(1) + elBeamSeparation:elBeamSeparation:elSweepRange(2);
            end
            azDir = repmat(azDir,1,beamsInEl);
            elDir = reshape(repmat(elDir,beamsInAz,1),[],1)';
            
        elseif numel(unique(azSweepRange)) == 2 && numel(unique(elSweepRange)) == 1 
            % Azimuth sweep only
            beamSeparation = (azSweepRange(2) - azSweepRange(1))/numBeams;
            azDir = azSweepRange(1) + beamSeparation:beamSeparation:azSweepRange(2);
            elDir = elSweepRange(1)*ones(1,numBeams);

        elseif numel(unique(azSweepRange)) == 1 && numel(unique(elSweepRange)) == 2 
            % Elevation sweep only
            beamSeparation = (elSweepRange(2) - elSweepRange(1))/numBeams;
            azDir = azSweepRange(1)*ones(1,numBeams);
            elDir = elSweepRange(1) + beamSeparation:beamSeparation:elSweepRange(2);

        else
            % No sweep
            azDir = azSweepRange(1)*ones(1,numBeams);
            elDir = elSweepRange(1)*ones(1,numBeams);
        end
        
        % Form azimuth and elevation angle pairs for all beams
        beamAng = [azDir;elDir];
    end
end