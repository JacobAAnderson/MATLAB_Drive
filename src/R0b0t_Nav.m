

classdef R0b0t_Nav
    
    properties
        Wp
        Done
        Goal
        Robot
        Type
        PF
        Cov
        
    end
    
    methods
        
        function obj = R0b0t_Nav(robot)
            
            obj.Wp.Waypoints = [];
            obj.Wp.indx = 1;
            obj.Done = false;
            obj.Goal = [];
            obj.Type = 'dr';
            
            switch lower(robot)
                case 'auv', obj.Robot = 'auv';
                case 'asv', obj.Robot = 'asv';
                otherwise, warning("%s is not the robot you are looking for", robot)
            end
        end
        
        
        function obj = Set_Nav(obj, type)
            
            switch lower(type)
                case 'dr',   obj.Type = 'dr';
                case 'slam', obj.Type = 'slam';
                case 'pf',   obj.Type = 'pf';
                otherwise, warning("%s is an invaled navigation type", type)
            end
            
        end
        
        
        % ____ Waypoint Navigation _________________________________________________
        function obj = Load_Waypoints(obj, file)
            
            [~,name,~] = fileparts(file);
            obj.Wp.file = name;
            obj.Wp.Waypoints = csvread(file);
            obj.Wp.indx = 1;
            obj.Done = false;
            
            if strcmpi(obj.Robot,'asv')
                obj.Wp.Waypoints(:,3) = 0;
            else
                obj.Wp.Waypoints(obj.Wp.Waypoints(:,3) > 0, 3) = -obj.Wp.Waypoints(obj.Wp.Waypoints(:,3) > 0, 3);
            end
            
        end
        
        
        function obj = WaypointFollowing(obj, X , r)
            
            obj.Goal = obj.Wp.Waypoints(obj.Wp.indx,:);
            
            if numel(obj.Goal) < 3 || strcmpi(obj.Robot,'asv') , obj.Goal(3) = 0; end % Make Sure thing are the right lengths
            if numel(X) < 3, X(3) = 0; end
            
            dist = obj.Goal(1:2)' - X(1:2);
            
            if sqrt(dist' * dist) < r, obj.Wp.indx = obj.Wp.indx +1; end
            
            if obj.Wp.indx > size(obj.Wp.Waypoints,1)
                obj.Wp.indx = 1;
                obj.Done = true;
            end
                        
        end
        
        
        function obj = Add_PartilceFilter(obj, pf), obj.PF = pf; end 
       
        
        function obj = InformativePathPlanning(obj, map, robotPath, start, maxDist, seeds)
            
            visitedLayer = zeros(size(map.Map));
            visitedLayer(isnan(map.Map)) = NaN;
            
            map = map.Add_Layer("Visited", visitedLayer);
            map = map.Set_LayerValue("Visited", robotPath(:,1:2), true);
%             map.PlotLayer("Visited");
            
            ev = Evolutionalry(map, start, 10, maxDist, seeds);
            obj.Wp.Waypoints = ev.InformativePath(2000, false);

%             ev = RRT(map, start, maxDist);
%             obj.Wp.Waypoints = ev.InformativePath(50, false);
            
            
            obj.Wp.indx = 2; % Index 1 is where the vehicle is
            
            if size(obj.Wp.Waypoints,2) < 3, obj.Wp.Waypoints(:,3) = 0; end
            
%             fig = map.PlotMap;
%             obj.PlotWaypoints(fig, 'path')
            
%             map.PlotLayer_Surf('Uncertainty');
            
        end
        
        
        %____ Plotting Functons ___________________________________________
        % Plot Waypoints
        function PlotWaypoints(obj,fig, type)
            
            figure(fig)
            hold on
            
            if nargin < 3 || strcmpi(type,'points')
                scatter3(obj.Wp.Waypoints(:,1),   obj.Wp.Waypoints(:,2),   obj.Wp.Waypoints(:,3),   '*g')
                scatter3(obj.Wp.Waypoints(1,1),   obj.Wp.Waypoints(1,2),   obj.Wp.Waypoints(1,3),   'dg')
                scatter3(obj.Wp.Waypoints(end,1), obj.Wp.Waypoints(end,2), obj.Wp.Waypoints(end,3), 'or')
                
            elseif strcmpi(type,'path')
                plot3(obj.Wp.Waypoints(:,1), obj.Wp.Waypoints(:,2), obj.Wp.Waypoints(:,3), 'g')
                scatter3(obj.Wp.Waypoints(1,1),   obj.Wp.Waypoints(1,2),   obj.Wp.Waypoints(1,3),   'dg')
                scatter3(obj.Wp.Waypoints(end,1), obj.Wp.Waypoints(end,2), obj.Wp.Waypoints(end,3), 'or')
                
            else
                warning("%s Is An Invalid Plotting type", type)
            end
            
            hold off
            drawnow
            
        end
        
    end
    
end