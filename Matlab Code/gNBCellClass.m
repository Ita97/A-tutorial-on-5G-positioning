%% 5G NR gNB cell Class
classdef gNBCellClass
   properties
        coordinates = [45.4787, 9.2330]      % [latitude, longitude]
        xyz_cartesian = [0 0 0]              % [x, y, z]
        orientation = [0 0]             % azimuth between -180 and 180 (deg)
                                        % elevation between -90 and 90 (deg)
        azimuthRangeAngle {mustBePositive} = 120    % (deg)
        elevationRangeAngle {mustBePositive} = 180  % (deg)
        pci = 0
        fc = 3.7e9            % cental frequency (Hz)
        los = false;
   end 
   properties(SetAccess=protected)
        arraySize = [8 8]     % 1, 2, 4, 8 or 16
        lambda 
        antennaArray 
   end
   methods
       function obj = gNBCellClass(coord, pci, orientation, arraySize, fc)

        obj.lambda = physconst("LightSpeed")/obj.fc;
        if nargin == 1
            obj.coordinates = coord;
        end
        if nargin == 2
            obj.coordinates = coord;
            obj.pci = pci;
        end
        if nargin == 3
            obj.coordinates = coord;
            obj.pci = pci;
            obj.orientation = orientation;
        end
        if nargin == 4
            obj.coordinates = coord;
            obj.pci = pci;
            obj.orientation = orientation;
            obj.arraySize = arraySize;
        end
        if nargin == 5
            obj.coordinates = coord;
            obj.pci = pci;
            obj.orientation = orientation;
            obj.fc = fc;
            obj.arraySize = arraySize;
        end
        
      end
      function r = getAzimuthLimits(obj)
          % return an array of two values
         ang = round(obj.azimuthRangeAngle/2);
         lim_sx = obj.orientation(1)-ang;
         if lim_sx < (-180)
             lim_sx = lim_sx + 360;
         end
         lim_dx = obj.orientation(1)+ang;
         if lim_dx > 180
             lim_dx = lim_dx - 360;
         end
         r = [lim_sx lim_dx];
      end
      function r = getElevationLimits(obj)
          % return an array of two values
         ang = round(obj.elevationRangeAngle/2);
         lim_sx = obj.orientation(2)-ang;
         if lim_sx < (-90)
             lim_sx = lim_sx + 180;
         end
         lim_dx = obj.orientation(2)+ang;
         if lim_dx > 90
             lim_dx = lim_dx - 180;
         end
         r = [lim_sx lim_dx];
      end
      function r = getCoordinates(obj, coord_idx)
          len = length(obj);
          if nargin > 1
              r = zeros(len,1);
              for i = 1: len
                  r(i)=obj(i).coordinates(coord_idx);
              end
          else
              r = zeros(len,2);
              for i = 1: len
                  r(i,:)=obj(i).coordinates;
              end
          end
      end
      function r = getXYZ(obj, axis)
          len = length(obj);
          if nargin > 1
              r = zeros(len,1);
              for i = 1: len
                  r(i)=obj(i).xyz_cartesian(axis);
              end
          else
              r = zeros(len,3);
              for i = 1: len
                  r(i,:)=obj(i).xyz_cartesian;
              end
          end
      end
      function obj = set.xyz_cartesian(obj, xyz_coord)
          len = length(obj);
          if len==size(xyz_coord,1)
              for i = 1: len
                  obj(i).xyz_cartesian=xyz_coord(i,:);
              end
          else
              error("Different size!")
          end
      end
      function r = getLatitude(obj)
          len = length(obj);
          r = zeros(len,1);
          for i = 1: len
              r(i)=obj(i).coordinates(1);
          end
          
      end
      function r = getLongitude(obj)
          len = length(obj);
          r = zeros(len,1);
          for i = 1: len
              r(i)=obj(i).coordinates(2);
          end
      end
      function r = getOrientation(obj)
          len = length(obj);
          r = zeros(len,2);
          for i = 1: len
              r(i,:)=obj(i).orientation;
          end
      end
      function obj = set.fc(obj, fc)
          len = length(obj);
          for i = 1: len
              obj(i).fc=fc;
              obj(i).lambda=physconst("LightSpeed")/fc;

%               obj(i).antennaArray = phased.NRRectangularPanelArray('ElementSet', ...
%                 {phased.NRAntennaElement},'Size',[obj(i).arraySize, 1, 1], ...
%                 'Spacing', [0.5*obj(i).lambda, 0.5*obj(i).lambda,1, 1]);
          end
      end
      function obj = set.arraySize(obj, array)
          len = length(obj);
          arraySz = zeros(1,2);
          if length(array)==1
              arraySz = [array array];
          end
          if length(array)==2
              arraySz=array;
          end
          if length(array)>2
              error("The size of array needs to be a scalar or an array 1x2")
          end
          for i = 1: len
              obj(i).arraySize=arraySz;
%               obj(i).antennaArray = phased.URA(obj(i).arraySize, 0.5*obj(i).lambda, 'Element', ...
%                 phased.IsotropicAntennaElement('BackBaffled',true));
              obj(i).antennaArray = phased.NRRectangularPanelArray('ElementSet', ...
                {phased.NRAntennaElement},'Size',[arraySz, 1, 1], ...
                'Spacing', [0.5*obj(i).lambda, 0.5*obj(i).lambda, 1, 1]);
          end
      end
      function r = isInAzimuthRange(obj, theta)
          o = obj.orientation(1);
          if o<0
              o = o + 360;
          end
          if theta<0
              theta = theta+360;
          end
          ang = round(obj.azimuthRangeAngle/2);
          r = (theta >= (o-ang)) && (theta <= (o+ang));
      end
     
   end
end