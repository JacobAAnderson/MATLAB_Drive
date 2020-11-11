% Geotiff Class
% Jacob Anderson
% 9/25/2019

classdef Geotiff
    properties
        Image
        Data
    end
    
    properties (Access = private)
        Data_utm
    end
    
    methods (Static)
        
        % Constructor
        function obj = Geotiff, warning on backtrace
            
            obj.Image = [];
            obj.Data = [];
            
        end
        
    end
    
    
    methods
        
        % Get Geotiff From file
        function obj = FromFile(obj, file)
            [obj.Image, obj.Data] = geotiffread(file);                      % Read Geo-tiff file
            obj = obj.Check_GeoData;                                        % Check that the geotiff is formated properly
        end
        
        
        % Get Geotiff From Google Earth
        function obj = FromGoogle(obj, lat, lon, zoom)
            
            [XX, YY, M, Mcolor] = get_google_map(lat, lon, 'Zoom', zoom); % Get Map from Google--> Max Zoom: 19
            
            geoImage = ind2rgb(M,Mcolor);                           % Convert into an RGB image
            obj.Image = flipud(geoImage);
            
            [~, ~, utmZone] = deg2utm(lat, lon);                    % get utm zone
            
            Xlim = [min(XX(:)), max(XX(:))];
            Ylim = [min(YY(:)), max(YY(:))];
            utmZone = [utmZone; utmZone];
            
            [latLim, lonLim] = utm2deg(Xlim', Ylim', utmZone);
            
            obj.Data = georasterref('RasterSize',     size(geoImage), ...
                'LatitudeLimits', [latLim(1), latLim(2)],   ...
                'LongitudeLimits',[lonLim(1), lonLim(2)]);
        end
        
        
        % Get Geotiff From Web Map Server
        function obj = FromWMS(obj, lat, lon)
            [obj.Image, obj.Data] = GetMapfromWMS( lon, lat );
        end
        
        
        % Display Geotiff
        function fig =  Show(obj)
            
            if nargout == 1
                fig = figure('name','Geotiff','numbertitle','off');
            end
            
            geoshow(obj.Image, obj.Data)
        end
        
        
        function points = GetPoints(obj)
            
            figure('name','Geotiff Get Points','numbertitle','off')
            obj.Show;
            
            [x,y] = getpts;
            
            points = [x,y];
            
        end
        
        
        % Save as .tiff
        function ExportTiff(obj, file)
            geotiffwrite(file, obj.Image, obj.Data);
        end
        
    end
    
    methods (Access = private)
        
        % Make Sure the Geo tiff is formated properly
        function obj = Check_GeoData(obj)
            
            obj.Image = im2uint8(obj.Image);                                             % Convert raster image to uint8 data type for faster graphic rendering
            
            if isempty(obj.Data)
                disp('Image is not a geo-referanced tiff')
                return
                
            elseif ~isa(obj.Data,'map.rasterref.GeographicCellsReference')               % Check that geoData is in the right format for 'geoshow'
                
                obj.Data = georasterref('RasterSize',size(obj.Image), ...                 % Reformat geoData for 'geoshow'
                    'RasterInterpretation', obj.Data.RasterInterpretation, ...
                    'LongitudeLimits',      obj.Data.XWorldLimits, ...
                    'LatitudeLimits',       obj.Data.YWorldLimits, ...
                    'ColumnsStartFrom',     obj.Data.ColumnsStartFrom,...
                    'RowsStartFrom',        obj.Data.RowsStartFrom );
            else
                disp('Geo Data checks out')
            end
            
        end
    end
    
end
