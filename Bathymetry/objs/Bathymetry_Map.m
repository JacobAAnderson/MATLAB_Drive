% Bathymetry Map CLAss
% Jacob Anderson
% 9/16/2019

classdef Bathymetry_Map
    
    properties
        Elevation
        Variance
        Easting
        Northing
        origin
        uperRight
        UTM_Zone
        Logs
        Header
        
    end
    
    properties (Access = private)
        FullRes_Elevation
        FullRes_Variance
        FullRes_Easting
        FullRes_Northing
        dims
    end
    
    methods
        
        % Constructor
        function obj = Bathymetry_Map(file, station_ID), warning on backtrace
            
            obj.FullRes_Elevation = [];
            obj.FullRes_Variance  = [];
            obj.FullRes_Easting   = [];
            obj.FullRes_Northing  = [];
            obj.Elevation         = [];
            obj.Variance          = [];
            obj.Easting           = [];
            obj.Northing          = [];
            obj.origin            = [];
            obj.uperRight         = [];
            obj.UTM_Zone  = '';
            obj.Logs      = '';
            obj.Header    = '';
            
            if nargin < 1, return, end                                      % If No data is provided, return empty bathymetry map
            
            emData = EM_Data(file);
            
            if nargin > 1
                emData = emData.TideCorrection(station_ID);                 % Get the Tide Level that correspondes to the when the data was colected
            end
            
            obj = obj.MakeMap(emData.RawData);
            
        end
        
        
        % Create the Map
        function obj = MakeMap(obj, Data, filterSpecs)
            
            if isstruct(Data)                                               % Get Data out of data structure
                vehicle    = cat(1, Data.vehicle);
                bathymetry = cat(1, Data.bathymetry);
                
                in = vehicle(:,3) < 0.1;                                        % Filter out data points with low GPS accuracey
                
                latitude  = bathymetry(in,1);
                longitude = bathymetry(in,2);
                data      = bathymetry(in,3);
                
                tide = cat(1, Data.tide);
                tide = tide(in);
                data = data - tide;
                
            else
                latitude  = Data(:,1);
                longitude = Data(:,2);
                data      = Data(:,3);
            end
            
            [x,y,utmZone] = deg2utm(latitude, longitude);
            
 
            
            if nargin < 3
                res       =   1.0;  % Grid cell width
                %                 baseFloor = -60;    % Floor of the map --> This will defalt to the minimum data value if the argument is not passed into the Bathymetry function.
                window    =  30;    % Window size of EKF
                sig       =   2;    % Sigma threshold for low pass filter
            else
                res       = filterSpecs(1);
%                 baseFloor = filterSpecs(2);
                window    = filterSpecs(3);
                sig       = filterSpecs(4);
            end
            
            % Make Bathymetry Map
            disp("Making Bathymetry Map")
            [Elivation, XX,  YY,  Var] = EKF_2D( x, y, data,'GridCell', res,'Window', window, 'sigma', sig); % ,'Floor', baseFloor );
            disp("Done")
            
            % Assign Data to struct
            obj.Header = {'   header: Explains the entries in this data structure'
                          'elevation: Water Depth [m]';
                          ' variance: Uncertanty in Water Depth [m]';
                          '  easting: UTM [m]';
                          ' northing: UTMs [m]'
                          '  utmZone: UTM zone'};
            
            obj.FullRes_Elevation = Elivation;
            obj.FullRes_Variance  = Var;
            obj.FullRes_Easting   = XX;
            obj.FullRes_Northing  = YY;
            
            obj.Elevation = Elivation;
            obj.Variance  = Var;
            obj.Easting   = XX;
            obj.Northing  = YY;
            obj.UTM_Zone  = utmZone(1,:);
            
            if isstruct(Data)
                obj.Logs = { Data(:).log }';
            else
                obj.Logs  = 'Lake Nigh tHorse Regions 3, 4, 5, 7, 8';
            end
            
            obj.origin    = [min(XX(:)), min(YY(:))];
            obj.uperRight = [max(XX(:)), max(XX(:))];
            obj.dims      = fliplr(size(XX));
            
        end
        
        
        function obj = MakeMap_GP(obj, data_lon_lat, parim_lon_lat)
            
            % data_lon_lat and parim_lon_lat are expected to be [lon, lat]
            
            [data_utm(:,1),  data_utm(:,2), utmZone]  = deg2utm(data_lon_lat(:,2),  data_lon_lat(:,1));
            [parim_utm(:,1), parim_utm(:,2), ~]       = deg2utm(parim_lon_lat(:,2), parim_lon_lat(:,1));
            
            
            if nargin == 3
                origin_  = min(parim_utm);
                topRight = max(parim_utm);
                X = [data_utm; parim_utm];
                Y = [data_lon_lat(:,3); parim_lon_lat(:,3)];
            else
                origin_  = min(data_utm);
                topRight = max(data_utm);
                X = data_utm(:,1:2);
                Y = data_lon_lat(:,3);
            end
            
            x = origin_(1): topRight(1);
            y = origin_(2): topRight(2);
            
            [XX, YY] = meshgrid(x,y);
            
            Xnew = [XX(:),YY(:)];
            
            
            disp("Making Bathymetry map with Gaussian Process Regrestion")
            
            gprMdl = fitrgp(X, Y);

            [ypred,ysd] = predict(gprMdl,Xnew); 
            
            disp('Done')
            
            Elevation_ = reshape(ypred, size(XX));
            Var        = reshape(ysd,   size(XX));
            
            obj.Header = {'   header: Explains the entries in this data structure'
                          'elevation: Water Depth [m]';
                          ' variance: Uncertanty in Water Depth [m]';
                          '  easting: UTM [m]';
                          ' northing: UTMs [m]'
                          '  utmZone: UTM zone'};
            
                      
            obj.FullRes_Elevation = Elevation_;
            obj.FullRes_Variance  = Var;
            obj.FullRes_Easting   = XX;
            obj.FullRes_Northing  = YY;
            
            obj.Elevation = Elevation_;
            obj.Variance  = Var;
            obj.Easting   = XX;
            obj.Northing  = YY;
            obj.UTM_Zone  = utmZone(1,:);
            
            obj.origin    = origin_(1:2);
            obj.uperRight = topRight(1:2);
            obj.dims      = fliplr(size(XX));
        end
        
        
        
        function obj = MapSmoothing(obj, fun)
           
            fprintf("\n\nSmoothing Bathymetry Map\n\n")
            
            if nargin < 2
                fun = @(x) nanmean(x);
            end
            
            map = obj.FullRes_Elevation;
            
            sz = size(map);
            
            sw = SlidingWindow([10,10], size(map));
            
            for ii = 1:size(map,1)
                for jj = 1:size(map,2)
                    
                    [rows, cols] = sw.D2([ii,jj]);
                    
                    if all(isnan(map(rows, cols))), continue, end
                    
                    ind = sub2ind(sz,rows,cols);
                    obj.FullRes_Elevation(ii,jj) = fun(map(ind));
                    
                end
            end
            
            obj.Elevation = obj.FullRes_Elevation;
            
        end
        
        
        % Get the Center of the map
        function center = Center(obj)
            
            xx = obj.Easting( ~isnan(obj.Elevation));
            yy = obj.Northing( ~isnan(obj.Elevation));
            
            center = [nanmean(xx), nanmean(yy)];
        end
        
        
        % Get Depth reading at a certain location
        function depth = Depth(obj, X)
            idx = state2index(obj, X(:,1:2));
            depth = obj.Elevation(idx);
        end
        
        
        % Povide a Lower Resolution Map
        function obj = ReduceResolution(obj, ratio)
            
            obj.Elevation = imresize( obj.FullRes_Elevation, ratio, 'method', 'nearest');
            obj.Variance  = imresize( obj.FullRes_Variance,  ratio, 'method', 'nearest');
            obj.Northing  = imresize( obj.FullRes_Northing,  ratio, 'method', 'nearest');
            obj.Easting   = imresize( obj.FullRes_Easting,   ratio, 'method', 'nearest');
            
            obj.origin    = [min(obj.Easting(:)), min(obj.Northing(:))];
            obj.uperRight = [max(obj.Easting(:)), max(obj.Northing(:))];
            obj.dims      = fliplr(size(obj.Easting));
            
        end
        
        
        % Export Batymetry Map as STL
        function Eport(obj, fileName)
            
            [~,~,ext] = fileparts(fileName);
            
            switch ext
                case '.stl'
                    
                    stlwrite(fileName, obj.Easting, obj.Northing, obj.Elevation);
                    
                case '.h5'
                    
                    h5create(fileName,'/Easting',size(obj.Easting))
                    h5write(fileName,'/Easting',obj.Easting)
                    
                    h5create(fileName,'/Northing',size(obj.Easting))
                    h5write(fileName,'/Northing',obj.Northing)
                    
                    h5create(fileName,'/Elevation',size(obj.Elevation))
                    h5write(fileName,'/Elevation',obj.Elevation)
                
                otherwise
                    warning("Bad Extendtion")
                 
            end
        end
        
        
        % Import Bathymetry from STL
        function obj = ImportBathymetry(obj, file)
            
            [x, y, z, c] = stlread(file);
            
            figure('NumberTitle','off','Name','Origonal Object')                        % Show what came in
            axis equal
            patch(x, y, z, c, 'FaceAlpha', 1)
            xlabel('X')
            ylabel('Y')
            zlabel('Z')
            view(45,30)
            
        end
        
        
        function mapp = AsMapp(obj)
            mapp = Mapp(obj.Easting, obj.Northing, obj.Elevation);
        end
            
            
        
        % Plotting Functions ----------------------------------------------
        % Plot the 3D Model
        function figOut = Plot_3DModel(obj, varargin)
            
            XX = obj.Easting;
            YY = obj.Northing;
            
            Bathymetry = obj.Elevation;
            
            fig = figure('Name','Bathymetry Model','NumberTitle','off');
            hold on
            
            %-- Create Surface from the Bathymetry Data --
            surface(XX, YY, -Bathymetry,'edgecolor', 'none');
            colormap(copper);
            caxis([(min(-Bathymetry(:)* 1.2)), 0])
            light('Position',[0 0 3000],'Style','infinite');
            lighting gouraud
            material dull
            
            %-- Overlay Contours --
            maxDepth = floor(min(-Bathymetry(:)));
            
            step = maxDepth / 10;
            
            
            v = 0:step:maxDepth;
            [c,h] = contour3(XX, YY,-Bathymetry,v, 'k');
            v = 0:-10:maxDepth;
            clabel(c,h,v,'FontSize',8)
            h = colorbar;
            
            title('Interpolated Bathymetry')
            xlabel('Easting [utm]')
            ylabel('Northing [utm]')
            zlabel('Depth [meters]')
            
            ylabel(h, 'Depth [meters]')
            
            if nargin > 1
                view(varargin{1}, varargin{2}) % Render image at a certain angle for viewing
            end
            
            hold off
            
            if narout == 1
                figOut = fig;
            end
            
        end
        
        
        % Plot the Bathymetry as a Contour Map
        function figOut = PlotMap(obj, geotiff)
            
            fig = figure('Name','Bathymetry Map', 'numbertitle','off');
            
            
            if nargin > 1                                                   % Plot Map with Geotiff if available
                
                geotiff.Show;
                
                % Convert to lat lon for the geotiff
                x = obj.Easting(1,:)';
                y = obj.Northing(:,1);
                
                zone = obj.UTM_Zone;
                
                [lat, ~]  = utm2deg( repmat(x(1), size(y)), y, repmat(zone, size(y)) );
                [~, lon]  = utm2deg( x, repmat(y(1), size(x)), repmat(zone, size(x)) );
                
                [LON,LAT] = meshgrid(lon, lat);
                
            else
                LON = obj.Easting;
                LAT = obj.Northing;
            end
            
            
            Elivation = obj.Elevation;
            hold on
            
            maxDepth = floor(min(-Elivation(:)));
            step = maxDepth / 10;
            
            v = 0:step:maxDepth;
            [c,h] = contourf(LON, LAT, -Elivation, v, 'k');
            v = 0:-10:maxDepth;
            clabel(c,h,v,'FontSize',8)
            h = colorbar;
            caxis([min(-Elivation(:))*1.2, 1])
            
            xtickangle(45)
            xlabel('Longitude [deg]')
            ylabel('Latitude [deg]')
            
            ylabel(h, 'Depth [Meters]')
            
            hold off
            
            if nargout == 1
                figOut = fig;
            end
        end
        
        
        
        % Plot the Bathymetry as a Contour Map
        function figOut = Plot_Map_UTM(obj, fig)
            
            if nargin > 2, figure(fig), hold on,
            else, fig = figure('Name','Bathymetry Map', 'numbertitle','off');
            end
            
            v = 0:-2.5:floor(min(-obj.Elevation(:)));
            
            [c,h] = contourf(obj.Easting, obj.Northing, -obj.Elevation, v, 'k');
            
            v = 0:-10:floor(min(-obj.Elevation(:)));
            clabel(c,h,v,'FontSize',8)
            
            h = colorbar;
            caxis([min(-obj.Elevation(:))-10, 1])
            
            xtickangle(45)
            xlabel('Easting [UTM]')
            ylabel('Northing [UTM]')
            
            ylabel(h, 'Depth [Meters]')
            
            hold off
            
            if nargout == 1
                figOut = fig;
            end
            
        end
        
        
        
        % Overlay a geotiff onto the 3D Model
        function figOut = Plot_3DModel_on_Geotiff(~, geotiff)
            
            fig = figure('Name', 'Surface with image projected on top of it');
            
            geotiff.Show;
            
%             [m,n,~] = size(geotiff.Image);
%             
%             x = linspace(geoData.LongitudeLimits(1), geoData.LongitudeLimits(2), n);
%             y = linspace(geoData.LatitudeLimits(1),  geoData.LatitudeLimits(2), m);
%             
%             [X,Y] = meshgrid(x,y);
%             
%             Z = X*0;
%             
%             [a,b] = size(Bathymetry);
%             
%             for ii = 1:(a*b)
%                 x = LON(ii);
%                 y = LAT(ii);
%                 z = Bathymetry(ii);
%                 
%                 if LON(ii) <
%                     
%                     row    = round( (LAT(ii) - geoData.LatitudeLimits(1))  *m/ (geoData.LatitudeLimits(2)-geoData.LatitudeLimits(1)));
%                     column = round( (LON(ii) - geoData.LongitudeLimits(1)) *n/ (geoData.LongitudeLimits(2)-geoData.LongitudeLimits(1)));
%                     
%                 end
%                 
%                 % Plot Bathymetry with geotiff
%                 
%                 I = rot90(geoImage);   % Rotate the image 180 degrees to display with the botom of the image forward
%                 I = rot90(I);
%                 
%                 surf(X, Y, Z, I,'EdgeColor','none'); % Plot surface
            
            
            if nargout == 1
                figOut = fig;
            end
            
        end
        
        
        % Save Object Function --> used by save()
        function s = saveobj(obj)
            s.Elevation = obj.FullRes_Elevation;
            s.Variance  = obj.FullRes_Variance;
            s.Easting   = obj.FullRes_Easting;
            s.Northing  = obj.FullRes_Northing;
            s.UTM_Zone  = obj.UTM_Zone;
            s.Logs      = obj.Logs;
            s.Header    = obj.Header;
        end
        
    end
    
    
    methods(Static)
        
        % Load object function --> used by load()
        function obj = loadobj(s)
            
            if isstruct(s)
                
                newObj =  Bathymetry_Map;
                
                newObj.FullRes_Elevation = s.Elevation;
                newObj.FullRes_Variance  = s.Variance;
                newObj.FullRes_Easting   = s.Easting;
                newObj.FullRes_Northing  = s.Northing;
                
                newObj.Elevation = s.Elevation;
                newObj.Variance  = s.Variance;
                newObj.Easting   = s.Easting;
                newObj.Northing  = s.Northing;
                newObj.UTM_Zone  = s.UTM_Zone;
                newObj.Logs      = s.Logs ;
                newObj.Header    = s.Header;
                
                newObj.origin    = [min(s.Easting(:)), min(s.Northing(:))];
                newObj.uperRight = [max(s.Easting(:)), max(s.Northing(:))];
                newObj.dims      = fliplr(size(s.Easting));
                
                obj = newObj;
            else
                obj = s;
            end
        end
        
    end
    
    
    methods (Access = private)
        
        % Get locaitons' indecies on map
        function idx = state2index(obj, X)
            diffs = (X - obj.origin) ./ (obj.uperRight - obj.origin);
            indicies = round( obj.dims .* diffs )+1;
            
            % Make sure indeicies do not exceed array Boundaries
            indicies(indicies(:,1)< 1, 1) = 1;
            indicies(indicies(:,2)< 1, 2) = 1;
            
            indicies(indicies(:,1)> obj.dims(1), 1) = obj.dims(1);
            indicies(indicies(:,2)> obj.dims(2), 2) = obj.dims(2);
            
            row = indicies(:,2);
            col = indicies(:,1);
            
            idx = sub2ind([obj.dims(2), obj.dims(1)], row, col);
        end
        
    end
end