function gNBPos_xyz = from_lat_long_to_xyz(bsSite, refsystemSite, varargin)
% Convert site objects coordinate in bsSite (lat, lon, high) according to
% refsystemSite site objects (lat, lon, high) in an xyz space.
    len = length(bsSite);
    gNBPos_xyz = zeros(len,3);
    % z
    % high = elevation(bsSite);
    % high_ref = elevation(ref_system);
    % gNBPos_xyz(3) = high - high_ref;
%     d = distance(ref_system,bsSite);
%     [az,el] = angle(ref_system,bsSite);
%     az = az * pi/180;
%     el = el * pi/180;
%     % x
%     gNBPos_xyz(1) = d * sin(az) * cos(el);
%     % y
%     gNBPos_xyz(2) = d * cos(az) * cos(el);
%     % z
%     gNBPos_xyz(3) = d * sin(el);


    % Validate sites
    validateattributes(bsSite,{'txsite','rxsite'},{'nonempty'}, ...
        'angle','',1);
    validateattributes(refsystemSite,{'txsite','rxsite'},{'nonempty'}, ...
        'angle','',2);
    % Process optional name/value pairs
    p = inputParser;
    addOptionalPath = nargin > 2 && mod(numel(varargin),2);
    travelPath = 'euclidean';
    if addOptionalPath
        % Validator function is necessary for inputParser to allow string
        % option instead of treating it like parameter name
        p.addOptional('Path',travelPath,@(x)ischar(x)||isstring(x));
    end
    p.addParameter('Map', []);
    p.addParameter('SourceAntennaSiteCoordinates', []);
    p.addParameter('TargetAntennaSiteCoordinates', []);
    p.parse(varargin{:});
    
    % Get usingCartesian from CoordinateSystem validation or from pre-validated 
    % AntennaSiteCoordinates
    if isempty(p.Results.SourceAntennaSiteCoordinates)
        usingCartesian = rfprop.internal.Validators.validateCoordinateSystem(bsSite, refsystemSite);
    else
        usingCartesian = strcmp(p.Results.SourceAntennaSiteCoordinates.CoordinateSystem,'cartesian');
    end
    
    if addOptionalPath
        try
            travelPath = validatestring(p.Results.Path, {'euclidean','greatcircle'}, ...
                'angle','',3);
        catch me
            % Check if option is a match for deprecated "geodesic" option
            try
                travelPath = validatestring(p.Results.Path, {'geodesic'}, ...
                    'angle','',3);
            catch
                rethrow(me)
            end
        end
    end
    
    % Get antenna site coordinates
    if usingCartesian
        map = rfprop.internal.Validators.validateCartesianMap(p);
        coordSys = 'cartesian';
    else
        map = rfprop.internal.Validators.validateMapTerrainSource(p, 'angle');
        coordSys = 'geographic';
    end
    rfprop.internal.Validators.validateMapCoordinateSystemMatch(map, coordSys);
    refCoords = rfprop.internal.Validators.validateAntennaSiteCoordinates(...
        p.Results.TargetAntennaSiteCoordinates, refsystemSite, map, 'angle');

    for i = 1:len
        bsCoords = rfprop.internal.Validators.validateAntennaSiteCoordinates(...
            p.Results.SourceAntennaSiteCoordinates, bsSite(i), map, 'angle');
        [x,y,z] = positionGeographic(bsCoords, refCoords);

        gNBPos_xyz(i,1) = x;
        gNBPos_xyz(i,2) = y;
        gNBPos_xyz(i,3) = z;
    end


end