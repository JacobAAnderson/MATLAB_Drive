

classdef RRT
    
    properties
        Map
        Nodes
        MaxDist
        Idx
%     end
%    
%     
%     properties (Access = private)
        dynamics
        noise
        derivatives
        limits
    end
    
    
    methods
        
        function obj = RRT(map, maxDist)
            
            obj.Map     = map;
            obj.MaxDist = maxDist;
            obj.Nodes   = struct('XYZ', [], 'Parent', [], 'Dist', [], 'Area', [], 'Info', []);
            obj.Idx     = 2;
            
            
            % Default Dynamics
            obj.dynamics.A = zeros(4);              % Vehicle Dynamics Matric
            obj.dynamics.A(1,3) = 1;
            obj.dynamics.A(2,4) = 1;
            
            obj.dynamics.B = zeros(4,2);            % Control input matrix
            obj.dynamics.B(3,1) = 1;
            obj.dynamics.B(4,2) = 1;
            
            % Limets for Kinodynamic Planning
            obj.limits.v = [-1,1];                  % Velocity Limits
            obj.limits.a = [-1,1];                  % Acceleration Limits
            obj.limits.u = [-1,1];                  % Control Input Limets
            
        end
        
        
        function obj = SetDynamics(obj, prop, func), obj.dynamics.(prop) = func; end        % Function for use to specify system functions
        
        function obj = SetDerivative(obj, prop, func), obj.derivatives.(prop) = func; end   % Set Derivatives
        
        function obj = SetLimets(obj, prop, limits), obj.limits.(prop) = limits; end
        
        function obj = SetNoise(obj, prop, dist, mu, sig)                   % Set Noise Distributions
            obj.noise.(prop) = makedist(dist,'mu',mu,'sigma',sig);
        end
        
        
        function path = Plan(obj, start, goal, epochs, plotIO)              % RRT path planning
            
            if numel(start) < 3, start(3) = 0; end
            if numel(goal)  < 3, goal(3) = 0;  end
            
            if nargin < 4, epochs = 1000; plotIO = flase;
            elseif nargin < 5, plotIO = false; end
            
            % Root Node
            obj.Nodes(1).XYZ    = start;
            obj.Nodes(1).Parent = 0;
            obj.Nodes(1).Dist   = [];
            obj.Nodes(1).Area   = [];
            obj.Nodes(1).Info   = [];
            
            
            % Build Tree
            for a = 1: epochs
                
                if mod(a, 10) == 0
                    rp = goal;
                else
                    rp = obj.Map.RandomPoints(1);                           % Generate Random Point
                    if numel(rp) < 3, rp(3) = 0; end
                end
                
                [nearID, nearXYZ] = obj.FindClossest(rp);                   % Find the Closest Node
                
                newXYZ = obj.Stear(rp, nearXYZ);                           % Stear towards the new node
                
                if ~obj.EvalPath(newXYZ, nearXYZ), continue, end            % Make Sure the Path is good
                
                obj.Nodes(obj.Idx).XYZ    = newXYZ;                         % Path is good, Add node
                obj.Nodes(obj.Idx).Parent = nearID;
                
                obj.Idx = obj.Idx + 1;
                
                if plotIO                                                   % Plot new data
                    plot(newXYZ(1), newXYZ(2), '*r')
                    line([nearXYZ(1), newXYZ(1)], [nearXYZ(2), newXYZ(2)],'Color','r');
                    drawnow
                end
                
                
                if vecnorm(newXYZ - goal) <= obj.MaxDist && obj.EvalPath(goal, newXYZ)  % Check for goal
                    obj.Nodes(obj.Idx).XYZ    = goal;
                    obj.Nodes(obj.Idx).Parent = obj.Idx - 1;
                    
                    if plotIO
                        line([newXYZ(1), goal(1)], [newXYZ(2), goal(2)],'Color','r');
                        drawnow
                    end
                    
                    break
                end
            end
            
            path = obj.PathToNode(obj.Idx);
            
        end
        
        
        function [path, control] = Plan_kd(obj, start, goal, epochs, plotIO)
            
            if numel(start) < 4, start(3:4) = 0; end
            if numel(goal)  < 4, goal(3:4)  = 0; end
            
            if nargin < 4, epochs = 1000; plotIO = flase;
            elseif nargin < 5, plotIO = false; end
            
            % Root Node
            obj.Nodes(1).XYZ    = start;
            obj.Nodes(1).Parent = 0;
            obj.Nodes(1).Dist   = [];
            obj.Nodes(1).Area   = [];
            obj.Nodes(1).Info   = [];
            
            for a = 1: epochs
                
                rp = zeros(4,1);
                
                if mod(a, 10) == 0
                    rp = goal;                                              % Sample the Goal
                else
                    rp(1:2) = obj.Map.RandomPoints(1);                                              % Generate Random Point       
                    if obj.Map.Collition(rp(1:2)), continue, end                                    % Make Sure the point is valad
                    rp(3:4) = rand(1,2) * (obj.limits.v(2) - obj.limits.v(2)) + obj.limits.v(1);    % Generate random Velocity
                end
                
                
            end
        
        end
    end
    
    methods (Access = private)
        
        function [minID, xyz] = FindClossest(obj, p)                        % Find the Closest Node
            
            test = cat(1, obj.Nodes.XYZ);
            dist = vecnorm(p - test, 2, 2);
            
            minID = find( dist == min(dist(:)));
            xyz = obj.Nodes(minID).XYZ;
        end
        
        
        function newXYZ = Stear(obj, p_Rand, p_Node)                        % Stear towards the sample point
            
            if vecnorm(p_Rand - p_Node) <= obj.MaxDist                      % If the length of the vetor does not exceed the max length then return that point as the new node
                newXYZ = p_Rand;
                return
            end
            
            v = p_Rand - p_Node;                                            % Vector from the nearest node to the random point
            n = vecnorm(v);                                                 % Length of vencor
            v = v/n;                                                        % Normalized Vector
            
            l = v * obj.MaxDist;                                            % Scale the vector to the desired length
            
            newXYZ = p_Node + l;                                            % Add the scaled vector to the location of the nearest node
            
        end
        
        
        function tf = EvalPath(obj, newXYZ, nearXYZ)                        % Check proposed path for obstacles
            
            V = newXYZ - nearXYZ;
            D = vecnorm(V);
            V = V/D;
            
            % Check for pathes leaving the map and / or in collitions
            ii = 0 : 0.1: D;
            test = nearXYZ + ii'.*V;
            
            if any(isnan(obj.Map.Measure(test))) || obj.Map.Collition(test)
                tf = false;
            else
                tf = true;
            end
            
        end
        
        
        function path = PathToNode(obj, id)                                 % Get path from goal to start
            
            path = zeros(100,3);
            idx = 1;
            
            while id ~= 0
                
                path(idx,:) = obj.Nodes(id).XYZ;
                id = obj.Nodes(id).Parent;
                idx = idx + 1;
                
            end
            
            path(idx : 100,:) = [];
            
            path = flipud(path);                                            % Flip path so that it is ordered start to goal
            
        end
        
    end
end




