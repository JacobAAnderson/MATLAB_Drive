


classdef Mapp
    
    properties
        Map             % Primary Map data
        Layers          % Additional data layers
        Easting         % X coordinates
        Northing        % Y Coordinates
        Area            % Totoal map area [ units are those of the Easting and northing fields]
        GridSize        % Size of a single grid eliment [x, y]
        GridArea        % Area of a simgle grid element
        Name            % Name of the Map
        Perim           % Perimiter of the Map
        Obstacles       % Environmental Obstacles
    end
    
    properties (Access = private)
        origin          % Origen of the Coordinate frame
        uperRight       % Top Roght corner of the coordinater frame
        dims            % Size of the maps grid [col "x", row "y"]
    end
    
    methods
        
        % Create Map
        function obj = Mapp(XX, YY, map, name), warning on backtrace
            
            obj.Map       = map;
            obj.Easting   = XX;
            obj.Northing  = YY;
            obj.origin    = [min(XX(:)), min(YY(:))];
            obj.uperRight = [max(XX(:)), max(YY(:))];
            obj.dims      = fliplr(size(XX));
            
            obj = obj.Get_Area;
            obj = obj.Get_GridSize;
            
            if nargin > 3, obj.Name = name;
            else, obj.Name = [];
            end
            
        end
        
        
        % Add a data layer to the map
        function obj = Add_Layer(obj, name, layer)
            
            if nargin < 3, layer = zeros(size(obj.Map)); end
            
            if any(fliplr(size(layer)) ~= obj.dims)
                warning("%s Layer sized Does not match map", name)
                return
            end
            
            obj.Layers.(name) = layer;
            
        end
        
        
        % Set the value of a layer
        function obj = Set_LayerValue(obj, type, X, value)
            idx = state2index(obj, X);
            obj.Layers.(type)(idx) = value;
        end
        
        
        % Access the a value on the map
        function value = Measure(obj, X, name)
            idx = state2index(obj, X);
            
            if nargin < 3, value = obj.Map(idx);
            else, value = obj.Layers.(name)(idx);
            end
            
        end
        
        
        % Indicate if a point lyes within the parimiter of the map
        function tf = InMap(obj, X), tf = inpolygon(X(:,1), X(:,2), obj.Perim(:,1), obj.Perim(:,2)); end
        
        
        % Get Info Layer names
        function list = InfoLayers(obj), list = fieldnames(obj.Layers); end
        
        
        % Generate Random Points on the Map
        function rp = RandomPoints(obj, num)
            
            rp = zeros(num, 2);
            
            ii = 1;
            while ii <= num
                
                p = rand(1,2) .* (obj.uperRight - obj.origin) + obj.origin;
                
                if isnan(obj.Measure(p)), continue, end
                
                rp(ii,:) = p;
                
                ii = ii+1;
            end
        end
        
        
        % ___ Obstacles and Collition Checking ____________________________
        function obj = Add_Obsticle(obj, poly)
            
            if ~isfield(obj.Obstacles, 'Poly')
                ind = 1;
            else
                ind = size(obj.Obstacles,2) + 1;
            end
            
            obj.Obstacles(ind).Poly = poly;
            
        end
        
        
        
        function tf = Collition(obj, x)
            
            if ~isfield(obj.Layers, 'obstacles') && ~isfield(obj.Obstacles, 'Poly')
                warning("Obstacles have not been added to this Map")
                tf = true;
                return
                
            elseif isfield(obj.Layers, 'obstacles') && ~isfield(obj.Obstacles, 'Poly')
                
                B = bwboundaries(obj.Layers.obstacles,'holes');
                
                for ii = 1: size(B,1)
                    poly = B{ii};
                    poly = [obj.Easting(1,poly(:,2))', obj.Northing(poly(:,1),1)];
                    poly(2:2:end,:) = [];
                    
                    obj.Obstacles(ii).Poly = poly;
                end
                
            end
            
            tf = false(size(obj.Obstacles,2),1);
            
            for ii = 1: size(obj.Obstacles,2)
                poly = obj.Obstacles(ii).Poly;
                tf(ii) = any( inpolygon(x(:,1), x(:,2), poly(:,1), poly(:,2)) );
            end
            
            tf = any(tf);
            
        end
        
        
        % Check for collitions
        function [h, s] = CollitionDist(obj, x, sig)
            
            if ~isfield(obj.Layers, 'obstacles') && ~isfield(obj.Obstacles, 'Poly')
                h = [];
                s  = [];
                warning("Obstacles have not been added to this Map")
                return
                
            elseif isfield(obj.Layers, 'obstacles') && ~isfield(obj.Obstacles, 'Poly')
                
                B = bwboundaries(obj.Layers.obstacles,'holes');
                
                for ii = 1: size(B,1)
                    poly = B{ii};
                    poly = [obj.Easting(1,poly(:,2))', obj.Northing(poly(:,1),1)];
                    poly(2:2:end,:) = [];
                    
                    obj.Obstacles(ii).Poly = poly;
                end
                
            end
            
            circ = error_ellips(x, sig);
            
            h = false;
            s = Inf;
            
            for ii = 1: size(obj.Obstacles,2)
                poly = obj.Obstacles(ii).Poly;
                
                for scale = 0.0: 0.25 : 10
                    tf =  hit(circ, poly, scale);
                    if tf
                        h = true;
                        s = min(s,scale);
                        break
                    end
                end
                
            end
            
        end
  
        
        % ____ Stat -- Indecies -- Neighbors ______________________________
        % Get locaitons' indecies on map
        function idx = state2index(obj, X)
            diffs = (X(:,1:2) - obj.origin) ./ (obj.uperRight - obj.origin);
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
        
        
        % Get locaitons' indecies on map
        function X = index2state(obj, ind)
            
            [row,col] = ind2sub([obj.dims(2), obj.dims(1)], ind);
            
            diffs = [col,row]./obj.dims;
            
            X = diffs .* (obj.uperRight - obj.origin) + obj.origin;
            
        end
        
        
        % Get neighbors
        function Xs = Neighbors(obj, X, eight_Connected)
            
            diffs = (X(:,1:2) - obj.origin) ./ (obj.uperRight - obj.origin);
            indicies = round( obj.dims .* diffs )+1;
            
            % Make sure indeicies do not exceed array Boundaries
            indicies(indicies(:,1)< 1, 1) = 1;
            indicies(indicies(:,2)< 1, 2) = 1;
            
            indicies(indicies(:,1)> obj.dims(1), 1) = obj.dims(1);
            indicies(indicies(:,2)> obj.dims(2), 2) = obj.dims(2);
            
            row = indicies(:,2);
            col = indicies(:,1);
            
            Xs = [row+1, col;
                row-1, col;
                row,   col-1;
                row,   col+1];
            
            if nargin > 2 && eight_Connected
                Xs = [Xs; row+1, col+1; row+1, col-1;
                    row-1, col+1; row-1, col-1];
            end
            
            idx = sub2ind([obj.dims(2), obj.dims(1)], Xs(:,1), Xs(:,2));
            
            Xs = index2state(obj, idx);
            
            tf = InMap(obj, Xs);
            
            Xs = Xs(tf,:);
            
        end
        
        
        % ___ Plotting ____________________________________________________
        % Plot the map
        function fig = PlotMap(obj, fig, alpha)
            
            ZZ = zeros(size(obj.Map));
            ZZ(isnan(obj.Map)) = NaN;
            
            im = 10 * obj.Map/nanmax(obj.Map(:));
            
            im(isnan(im)) = 0;
            
            im = round(im);
            
            map = parula(10);
            rgbImage = ind2rgb(im, map);
            
            % Create Figure
            if nargin > 1, figure(fig); hold on                             % Use provided figure handel
            elseif isempty(obj.Name), fig = figure('name', 'Mapp');         % Create new figure, no name
            else, fig = figure('name', sprintf("Map of %s", obj.Name));     % Create new figure, map has a name
            end
            
            fig.NextPlot = 'replacechildren';
            
            m = surf(obj.Easting, obj.Northing, ZZ, rgbImage, 'EdgeColor', 'none');
            view(0,90)
            
            if nargin > 2
                m.FaceAlpha = alpha;
                m.EdgeAlpha = alpha;
            end
            
            if isfield(obj.Obstacles, 'Poly')
                hold on
                
                for ii = 1: size(obj.Obstacles,2)
                    poly = obj.Obstacles(ii).Poly;
                    plot(poly(:,1), poly(:,2))
                end
                
            end
            
            hold off
            
            if nargin < 2
                title(sprintf("Map of %s", obj.Name));
                xlabel('Easing')
                ylabel('Northing')
            end
            
        end
        
        
        % Plot Map Parimiter
        function fig = PlotPerimiter(obj, fig, trans)
            
            % Create Figure
            if nargin > 1, figure(fig); hold on                             % Use provided figure handel
            elseif isempty(obj.Name), fig = figure('name', 'Mapp');         % Create new figure, no name
            else, fig = figure('name', sprintf("Map of %s", obj.Name));       % Create new figure, map has a name
            end
            
            fig.NextPlot = 'replacechildren';
            
            p = plot(obj.Perim(:,1), obj.Perim(:,2));
            
            if nargin > 2, alpha(p, trans); end
            
            hold off
            
            if nargin < 2
                title(sprintf("Map of %s", obj.Name));
                xlabel('Easing')
                ylabel('Northing')
            end
            
        end
        
        
        % Plot information layer
        function fig = PlotLayer(obj, type, fig, alpha)
            
            if ~isfield(obj.Layers,type)
                fig = [];
                warning("%s is not an information layer in this Map", type)
                return
            end
            
            ZZ = zeros(size(obj.Map));
            ZZ(isnan(obj.Map)) = NaN;
            
            im = 10 * obj.Layers.(type) / nanmax(obj.Layers.(type)(:));
            
            im(isnan(im)) = 0;
            
            im = round(im);
            
            map = parula(10);
            rgbImage = ind2rgb(im, map);
            
            % Create Figure
            if nargin > 2, figure(fig); hold on                                         % Use provided figure handel
            elseif isempty(obj.Name), fig = figure('name', sprintf("Map of %s", type)); % Create new figure, map has a name
            else, fig = figure('name', sprintf("Map of %s: %s", obj.Name, type));       % Create new figure, map has a name
            end
            
            fig.NextPlot = 'replacechildren';
            
            m = surf(obj.Easting, obj.Northing, ZZ, rgbImage, 'EdgeColor', 'none');
            view(0,90)
            
            if nargin > 3
                m.FaceAlpha = alpha;
                m.EdgeAlpha = alpha;
            end
            
            hold off
            
            if nargin < 3
                if isempty(obj.Name), title(sprintf("Map of %s", obj.Name));
                else, title(sprintf("Map of %s: %s", obj.Name, type));
                end
                
                xlabel('Easting')
                ylabel('Northing')
            end
            
        end
        
        
        % Plot information layer
        function fig = PlotLayer_Surf_WithMap(obj, type, fig, alpha)
            
            if ~isfield(obj.Layers,type)
                fig = [];
                warning("%s is not an information layer in this Map", type)
                return
            end
            
            im = 10 * obj.Map/nanmax(obj.Map(:));
            
            im(isnan(im)) = 0;
            
            im = round(im);
            
            map = parula(10);
            rgbImage = ind2rgb(im, map);
            
            % Create Figure
            if nargin > 2, figure(fig); hold on                                         % Use provided figure handel
            elseif isempty(obj.Name), fig = figure('name', sprintf("Map of %s", type)); % Create new figure, map has a name
            else, fig = figure('name', sprintf("Map of %s: %s", obj.Name, type));       % Create new figure, map has a name
            end
            
            m = surf(obj.Easting, obj.Northing, obj.Layers.(type), rgbImage, 'EdgeColor', 'none');
            hold on
            
            levels = linspace(nanmin(obj.Layers.(type)(:)), nanmax(obj.Layers.(type)(:)), 10);
            contour3(obj.Easting, obj.Northing, obj.Layers.(type), levels, 'r');
            
            view(0,70)
            
            if nargin > 3
                m.FaceAlpha = alpha;
                m.EdgeAlpha = alpha;
            end
            
            hold off
            
            if nargin < 3
                if isempty(obj.Name), title(sprintf("Map of %s", obj.Name));
                else, title(sprintf("Map of %s: %s", obj.Name, type));
                end
                
                xlabel('Easting')
                ylabel('Northing')
            end
            
        end
        
        
        % Plot information layer
        function fig = PlotLayer_Surf(obj, type, fig, alpha)
            
            if ~isfield(obj.Layers,type)
                fig = [];
                warning("%s is not an information layer in this Map", type)
                return
            end
            
            % Create Figure
            if nargin > 2, figure(fig); hold on                                         % Use provided figure handel
            elseif isempty(obj.Name), fig = figure('name', sprintf("Map of %s", type)); % Create new figure, map has a name
            else, fig = figure('name', sprintf("Map of %s: %s", obj.Name, type));       % Create new figure, map has a name
            end
            
            m = surf(obj.Easting, obj.Northing, obj.Layers.(type), 'EdgeColor', 'none');
            hold on
            
            levels = linspace(nanmin(obj.Layers.(type)(:)), nanmax(obj.Layers.(type)(:)), 10);
            contour3(obj.Easting, obj.Northing, obj.Layers.(type), levels, 'k');
            
            view(0,70)
            
            if nargin > 3
                m.FaceAlpha = alpha;
                m.EdgeAlpha = alpha;
            end
            
            hold off
            
            if nargin < 3
                if isempty(obj.Name), title(sprintf("Map of %s", obj.Name));
                else, title(sprintf("Map of %s: %s", obj.Name, type));
                end
                
                xlabel('Easting')
                ylabel('Northing')
            end
            
        end
        
        
        % Plot information layer
        function fig = PlotLayer_Contour(obj, type, fig)
            
            if ~isfield(obj.Layers,type)
                fig = [];
                warning("%s is not an information layer in this Map", type)
                return
            end
            
            % Create Figure
            if nargin > 2, figure(fig); hold on                                         % Use provided figure handel
            elseif isempty(obj.Name), fig = figure('name', sprintf("Map of %s", type)); % Create new figure, map has a name
            else, fig = figure('name', sprintf("Map of %s: %s", obj.Name, type));       % Create new figure, map has a name
            end
            
            levels = linspace(nanmin(obj.Layers.(type)(:)), nanmax(obj.Layers.(type)(:)), 10);
            contourf(obj.Easting, obj.Northing, obj.Layers.(type), levels, 'k');
            
            hold off
            
            if nargin < 3
                if isempty(obj.Name), title(sprintf("Map of %s", obj.Name));
                else, title(sprintf("Map of %s: %s", obj.Name, type));
                end
                
                xlabel('Easting')
                ylabel('Northing')
            end
            
        end
        
        
        % Plot information layer
        function fig = PlotMap_Contour(obj,fig)
            
            % Create Figure
            if nargin > 1, figure(fig); hold on                                         % Use provided figure handel
            else, fig = figure('name', sprintf("Map of %s", obj.Name));   % Create new figure, map has a name
            end
            
            levels = linspace(nanmin(obj.Map(:)), nanmax(obj.Map(:)), 10);
            contourf(obj.Easting, obj.Northing, obj.Map, levels, 'k');
            colorbar
            
            hold off
            if nargin < 2
                xlabel('Easting')
                ylabel('Northing')
                title(sprintf("Map of %s", obj.Name))
            end
            
        end
        
        
    end
    
    
    
    methods (Access = private)
        
        % Calculate the size of an individual grid elemant
        function obj = Get_GridSize(obj)
            
            diff = obj.uperRight - obj.origin;
            
            obj.GridSize = diff ./ fliplr(size(obj.Map));
            
            obj.GridArea = obj.GridSize(1) * obj.GridSize(2);
            
        end

        
        % Calculate the ara of an individual grid element
        function obj = Get_Area(obj)
            
            BW = obj.Map;
            BW(isnan(BW)) = 0;
            BW = logical(BW);
            
            B = bwboundaries(BW);
            B = B{1};
            
            idx = sub2ind([obj.dims(2), obj.dims(1)], B(:,1), B(:,2));
            
            perrim = [obj.Easting(idx), obj.Northing(idx)];
            
            obj.Area = polyarea(perrim(:,1), perrim(:,2));
            obj.Perim = perrim;
            
            %             perim = false(size(BW));                                      % Occupancy Grid ??
            %             perim(idx) = true;
            %
            %             obj = obj.Add_Layer('Perimiter', perim);
            
        end
        
        
    end
    
end




function tf = hit(circ, poly, scale)

x = mean(circ, 1);

S = [scale, 0, 0; 0, scale, 0; 0, 0, 1];

for ii = 1: size(circ,1)
    tf = [circ(ii,:),1] * S;
    circ(ii,:) = tf(1:2);
end

circ = circ - mean(circ, 1) + x;

tf = any( inpolygon(circ(:,1), circ(:,2), poly(:,1), poly(:,2)) );

end


function circ = error_ellips(x, sig)

circ = [cos(0:0.1:2*pi)', sin(0:0.1:2*pi)'];

R = [chol(sig), [0;0]; x, 1];

for ii = 1: size(circ,1)
    tf = [circ(ii,:),1] * R;
    circ(ii,:) = tf(1:2);
end

end
