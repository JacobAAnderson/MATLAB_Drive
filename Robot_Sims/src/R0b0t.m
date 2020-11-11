% Robot Class for Sims
% Jake Anderson
% 11/6/2019


classdef R0b0t
    
    properties
        Name
        Type
        State
        Paths
        Sensors
        Reward
        Comms
        Models
        Nav
        Goal
        RecivedData
        EKF
        Planning
%     end
%     
%     properties (Access = private)
        Map
        Noise                                                               % Noise Models
        Dynamics
        Limits
    end
    
    
    methods     
%% Set Up

        % Constructor
        function obj = R0b0t(type, name)
            
            switch type
                
                case 'ASV'
                    obj.Type = 'asv';
                    
                case 'AUV'
                    obj.Type = 'auv';
                    
                case 'Glider'
                    obj.Type = 'glider';
                    
                otherwise
                    warning('Un-Available Robot Class')
            end
            
            % Models
            obj.Nav = R0b0t_Nav(type);
            
            obj.Dynamics.f = @(x,u) x + u;                                  % Default vehicle dynamics
            
            obj.Noise.Speed   = makedist('Normal','mu',0,'sigma',0.25);     % Default vehilce transition model
            obj.Noise.Heading = makedist('Normal','mu',0,'sigma',5*pi/180);        % Default vehilce transition model
%             obj.Noise.Pitch   = makedist('Normal','mu',0,'sigma',3);
            
            obj.State.XYZ.DR  = [0, 0, 0, 0]';                              % Deadreckoning Vehicle State [ x, y, heading, speed]
            obj.State.XYZ.GT  = [0, 0, 0, 0]';                              % Ground Trouth Vehicle State [ x, y, heading, speed]
            obj.State.Cov = eye(4);
            
            obj.Limits.maxA = 2;                                            % Max Acceleration
            obj.Limits.maxV = 5;                                            % Max Velocity
            obj.Limits.maxTurn = 30 * pi/180;                               % Max Turn Angle 
            
            obj.Paths.GT  = zeros(1000,4);
            obj.Paths.DR  = zeros(1000,4);
            obj.Paths.PF  = zeros(1000,4);
            obj.Paths.PF_indx  = 1;
            obj.Paths.ind = 1;
            
            obj.Reward = 0;
            
            if nargin > 1, obj.Name = name; end
            
        end
        
        
        % Add Global Coordinates
        function obj = Add_Map(obj, map), obj.Map = map; end
        
        
        % Place Robt at certain place
        function obj = Set_VehicleState(obj, X)
            
            if numel(X) < 3, X = [X, 0, 0]; end
            
            obj.State.XYZ.DR = X';                                          % Set State
            obj.State.XYZ.GT = X';
            
            obj.Paths.DR(obj.Paths.ind,:) = obj.State.XYZ.DR;               % Show change in the path
            obj.Paths.GT(obj.Paths.ind,:) = obj.State.XYZ.GT;
        end
        
        
        % Span Robot in random location
        function obj = Spawn(obj, xlim, ylim)
            
            x = rand(4) .* [(max(xlim) - min(xlim)) + min(xlim);
                            (max(ylim) - min(ylim)) + min(ylim);
                                                            360;
                                                              3];
            obj.State.XYZ.DR = x;
            obj.State.XYZ.GT = x;
            
        end
        
        
        % Function for use to specify system functions
        function obj = SetDynamics(obj, prop, func), obj.Dynamics.(prop) = func; end 
        
        
        % Set vehicle transition model
        function obj = SetUncertainty(obj, type, dist, mu, sig)
            obj.Noise.(type) = makedist(dist,'mu',mu,'sigma',sig);
        end
        

        
%% Kinematicks
       
        % Move the vehicel towards a goal
        function obj = Move(obj, goal, speed, dt)
            
            if numel(goal) < 4, goal = [goal(1:2), 0, 0]; end
            
            dist = goal' - obj.State.XYZ.DR;                                 % Distance between robot and goal in cartesian spcae
            
            if sqrt(dist(1:2)' * dist(1:2)) < 5 && dist(3) * dt > 4 
                fprintf('\n%s Slowing Down!!\n', obj.Name)
                speed = 1; 
            end     % If we're there, stop
            
            bearing = atan2(dist(2), dist(1));                               % bearing to the goal
            
            % Signed shortes angle between the bearing and heading
            phi = (bearing - obj.State.XYZ.DR(3));                          
            
            if abs(phi) > 180
                phi = mod( mod((phi),2*pi)+540 * pi/180, pi ) - pi;
            end
            
            if abs(phi) > obj.Limits.maxTurn 
                phi = obj.Limits.maxTurn * sign(phi);
            end
            
            % Control inputs
            a = min((speed - obj.State.XYZ.DR(4)) / dt, obj.Limits.maxA);
            u =  0.25 .* [a; phi];
            
            % ___ State update ____________________________________________
            %     --> what the robot thinks is happening
            
            if isa(obj.EKF, 'ExtendedKalmanFilter') && isa( obj.Nav.PF, 'ParticleFilter')
                
                z = obj.Nav.PF.X;
                R = obj.Nav.PF.Cov;
                
                z(3) = CompassAngle(z(3), 'deg', 'rad');                    % Conver from Compas angle in degrres to Z up angle in radians
                
                x = obj.State.XYZ.DR;
                
                if abs(x(3) - z(3)) > pi
                    [x(3), z(3)] = MinAngle(x(3), z(3));
                end
                
                [X2, sig2] = obj.EKF.Update(x, obj.State.Cov, u, z, dt, R);
                obj.State.XYZ.DR = X2;
                obj.State.Cov = sig2;
                
                
                % Show Vehicle Covariance
                %                 fprintf('\n\n%s\n',obj.Name)
                %                 disp(sig2)
                
                
            else
                obj.State.XYZ.DR = obj.Dynamics.f(obj.State.XYZ.DR, u, dt, [0;0]);
            end
            
            %     --> What's Really Happening 
            n = [ random(obj.Noise.Heading, 1); random(obj.Noise.Speed,1)];
            
            switch obj.Type
                
                case 'auv'
                    obj.State.XYZ.GT = obj.Dynamics.f(obj.State.XYZ.GT, u, dt, n);
                    
                case 'glider'
                    obj.State.XYZ.GT = obj.State.XYZ.GT + [dX, dY, dZ];
                    obj.State.XYZ.GT = obj.Dynamics.f(obj.State.XYZ.GT, u, dt, n);
                    
                case 'asv'
                    obj.State.XYZ.GT = obj.Dynamics.f(obj.State.XYZ.DR, u, dt, n);
                    
                otherwise
                    warning('This is not the vehicle you are looking for')
            end
            
            % ____ Add Movement to path ___________________________________
            obj.Paths.DR(obj.Paths.ind,:) = obj.State.XYZ.DR;
            obj.Paths.GT(obj.Paths.ind,:) = obj.State.XYZ.GT;
            obj.Paths.ind = obj.Paths.ind + 1;
            
        end
        
        
        % has the Robot Made it to the goal
        function tr = Goal_Achived(obj, r)
            
            if numel(obj.Goal) < 4, obj.Goal = [obj.Goal(1:2), 0, 0]; end
            
            dist = obj.Goal' - obj.State.XYZ.DR;
            
            if sqrt(dist' * dist) < r, tr = true;
            else, tr = false;
            end
            
        end
        
        

%% Path Planning
        function obj = Step(obj, sim_time)
            
            if mod(sim_time, 60) == 0 || obj.Nav.Done
                fprintf("Replanning: %f\n", sim_time)
                obj = obj.MakeWaterPropertyMap("Salt", false);
                obj = obj.PathPlanning;
                obj.Nav.Done = false;
            end
            
        end
        
        
        function obj = PathPlanning(obj, seeds)
            
            map = obj.Models.("Salt").Map;
            
            path = obj.Paths.DR(:,1:2);
            path(path(:,1) == 0, :) = [];
            
            maxDist = 4*60;
            
            if nargin > 1
                obj.Nav = obj.Nav.InformativePathPlanning(map, path, obj.State.XYZ.DR, maxDist, seeds);
            else
                obj.Nav = obj.Nav.InformativePathPlanning(map, path, obj.State.XYZ.DR, maxDist);
            end
            
        end
           
        
        
%% Navigation Interface

        function obj = Navigate(obj,speed, dt, r)
            
            obj.Nav = obj.Nav.WaypointFollowing(obj.State.XYZ.DR(1:2), r);
            
            obj = obj.Move(obj.Nav.Goal, speed, dt);
            
        end

        
        function obj = UpdatePF(obj, dt, add_sensors)
            
            alt = obj.Sensors.Depth.value( obj.Sensors.Depth.iter );
            
            % heading = obj.State.XYZ.DR(3)*180/pi;
            heading = obj.Sensors.Compass.value(obj.Sensors.Compass.iter) * 180/pi;
            
            if heading < 90, heading = 90 - heading;
            else, heading = 360 - (heading - 90);
            end
            
            if nargin < 3
                obj.Nav.PF = obj.Nav.PF.Update(obj.State.XYZ.DR(4), dt, heading, alt);
            else
                obj.Nav.PF = obj.Nav.PF.Update_Acoms(obj.State.XYZ.DR(4), dt, heading, alt, add_sensors);
            end
            
            obj.Paths.PF(obj.Paths.PF_indx,:) = obj.Nav.PF.X;
            
            obj.Nav.Cov(obj.Paths.PF_indx).PF = obj.Nav.PF.Cov;
            
            obj.Paths.PF_indx = obj.Paths.PF_indx + 1;
        end
        
       
        
%% Data Collection and Processing

        % Create a sensor
        function obj = Add_Sensor(obj, type, noise_dist, mu, sig)
            
            if nargin < 5
                noise_dist = 'Normal';
                mu = 0;
                sig = 0;
            end
            
            obj.Sensors.(type).value = [];
            obj.Sensors.(type).xyz   = obj.State.XYZ.DR';
            obj.Sensors.(type).time  = [];
            obj.Sensors.(type).noise = makedist(noise_dist,'mu',mu,'sigma',sig);
            obj.Sensors.(type).iter  = 0;
            
        
        end
        
        
        % Add Sensor measurments
        function obj = Add_SensorMeasurment(obj, type, value, time)
            
            ii = obj.Sensors.(type).iter + 1;
            
            obj.Sensors.(type).value(ii,:) = value + random(obj.Sensors.(type).noise, 1,numel(value));
            obj.Sensors.(type).xyz(ii,:) = obj.State.XYZ.DR';
            obj.Sensors.(type).time(ii)  = time;
            obj.Sensors.(type).iter      = ii;
        end
        
        
             
%% Communications

        % Add comunication mode
        function obj = Add_Comms(obj, type)
            obj.Comms.(type).message = {};
            obj.Comms.(type).time = [];
            obj.Comms.(type).iter = 0;        
        end

        
        % Recive Communications
        function obj = Recive_Comms(obj, type, message, time)
            ii = obj.Comms.(type).iter + 1;
            obj.Comms.(type).message(ii,:) = message;
            obj.Comms.(type).time = time;
            obj.Comms.(type).iter = ii;
        end
        
    
        
%% Water Science Stuff 

        % Create Water Property Map
        function obj = MakeWaterPropertyMap(obj, type, plotIO)
            
            X = cat(1,obj.Sensors.(type).xyz);
            y = [obj.Sensors.(type).value];
            
            obj.Models.(type) = WaterFeature(obj.Map.Easting, obj.Map.Northing, obj.Map.Map, type);
            
            if strcmp(obj.Type,'asv')                                       % Save some space and time
                X = X(:,1:2);
                obj.Models.(type) = obj.Models.(type).MakeFromData(X, y);
                
                if nargin > 2 && plotIO, obj.Models.(type).PlotMap; end
                
            else
                obj.Models.(type) = obj.Models.(type). Make3DModelFromData(X, y, obj.Map.Map);
            end
            
        end
        
        
        % Predict Features
        function obj = PredictFeatures(obj, type, plotIO)
            
            X = cat(1,obj.Sensors.(type).xyz);
            y = [obj.Sensors.(type).value];
            
            obj.Models.(type) = WaterFeature(obj.Map.Easting, obj.Map.Northing, obj.Map.Map);
            
            if strcmp(obj.Type,'asv')                                       % Save some space and time
                X = X(:,1:2);
                obj.Models.(type) = obj.Models.(type).MakeFromData(X, y);
                obj.Models.(type) = obj.Models.(type).Find_Gradiant;
                obj.Models.(type) = obj.Models.(type).Fit_Model;
                
                if nargin > 2 && plotIO
                    obj.Models.(type).PlotGP_Prediction;
                    fig = obj.Models.(type).PlotMap;
                    obj.Models.(type).Plot_MaxMin([0,90], fig);
                    obj.Models.(type).PlotGrad;
                    obj.Models.(type).PlotModel;
                end
                
            else
%                 obj.Models.(type) = obj.Models.(type). Make3DModelFromData(X, y, obj.Map.Map);
                
                X = X(:,1:2);
                obj.Models.(type) = obj.Models.(type).MakeFromData(X, y);
                
                if nargin > 2 && plotIO
                    obj.Models.(type).PlotMap;
                    
                end
            end
            
        end
        
        
        % Transmit Predicted Water Features
        function objects = SendWaterPropertyMap(obj, type)
            
            objects = obj.Models.(type).Features.property;
%             fprintf("\n\n%s Sending %s Map\n\n", obj.Name, type)
            
        end
        
        
        % Recive Water Features
        function obj = ReciveWaterPropertyMap(obj, type, objects)
            
            obj.RecivedData.(type) = objects;
%             fprintf("\n\n%s Reciving %s Map\n\n", obj.Name, type)
            
            obj.Models.(type) = WaterFeature(obj.Map.Easting, obj.Map.Northing, obj.Map.Map, type);
            
            for ii = 1: size(obj.RecivedData.(type),2)
                
                mu    = obj.RecivedData.(type)(ii).mu; 
                sigma = obj.RecivedData.(type)(ii).sigma;
                rho   = obj.RecivedData.(type)(ii).r;
                s     = obj.RecivedData.(type)(ii).s;
                
                obj.Models.(type) = obj.Models.(type).Add_Gaussian(mu, sigma, rho, s);
                
                
            end
            
            obj.Models.(type) = obj.Models.(type).MakeMap_world;
%             obj.Models.(type).Map.PlotMap_Contour
            
            obj = obj.PathPlanning(objects);
            
        end
        

        
%% Stats on data

        % Compute Error between DR location and Ground truth
        function stats = PathStates(obj, plotIO)
            
            err = obj.Paths.DR(:,1:2) - obj.Paths.GT(:,1:2);
            stats.err = vecnorm(err,2,2);
            
            if plotIO
                figure('name', 'Path States', 'numbertitle', 'off')
                plot(stats.err)
                title(sprintf('\n\n%s Path RMSE\n\n', obj.Name))
                xlabel('time [s]')
                ylabel('RMSE [m]')
            end
                
            
        end
        
        
%% Plotting 
               
        % Plot the location of the robot
        function PlotLocation(obj, fig)
            
            figure(fig)
            hold on
            switch obj.Type
                
                case 'auv'
                    
                    s1 = scatter(obj.State.XYZ.DR(1), obj.State.XYZ.DR(2), '.c', ...
                        'MarkerEdgeAlpha', 0.6, ...
                        'MarkerFaceAlpha', 0.6);
                    s2 = scatter(obj.State.XYZ.GT(1), obj.State.XYZ.GT(2), '.r');
                    legend([s1,s2], 'DR', 'GT')
                    
                case 'asv'
                    scatter(obj.State.XYZ.GT(1), obj.State.XYZ.GT(2), '.g')
%                     scatter3(obj.State.XYZ.DR(1), obj.State.XYZ.DR(2), obj.State.XYZ.DR(3), '.c', ...
%                         'MarkerEdgeAlpha', 0.6, ...
%                         'MarkerFaceAlpha', 0.6)
                    
                otherwise
                    warning("We Can Not Disclose the locaiton of this device")
            end
            
            hold off
            drawnow
        end
        
        
        % Plot GT and DR paths
        function PlotPaths(obj, fig)
            
            i = obj.Paths.ind - 1;
            
            if nargin > 1, figure(fig); hold on
            else, figure;
            end
            
            switch obj.Nav.Type
                
                case 'dr' % Plot with Partilce filter
                    p1 = plot3(obj.Paths.GT(1:i,1), obj.Paths.GT(1:i,2), ones(i,1), 'g'); % obj.Paths.GT(1:i,3), 'c')
                    hold on
                    p2 = plot3(obj.Paths.DR(1:i,1), obj.Paths.DR(1:i,2), ones(i,1), 'r'); %obj.Paths.DR(1:i,3), 'r')
%                     p3 = plot3(obj.Paths.PF(1:i,1), obj.Paths.PF(1:i,2), ones(i,1), 'b'); %obj.Paths.DR(1:i,3), 'b')
                    legend([p1, p2], 'GT', 'DR')
                    
                case 'pf' % Plot with Partilce filter
                    p1 = plot3(obj.Paths.GT(1:i,1), obj.Paths.GT(1:i,2), ones(i,1), 'g'); % obj.Paths.GT(1:i,3), 'c')
                    hold on
%                     p2 = plot3(obj.Paths.DR(1:i,1), obj.Paths.DR(1:i,2), ones(i,1), 'r'); %obj.Paths.DR(1:i,3), 'r')
                    p3 = plot3(obj.Paths.PF(1:i,1), obj.Paths.PF(1:i,2), ones(i,1), 'b'); %obj.Paths.DR(1:i,3), 'b')
%                     legend([p1, p2, p3], 'GT', 'DR', 'PF')
                    legend([p1, p3], 'GT', 'TBN')
                    
                otherwise
                    warning("We Can Not Disclose the locaiton of this device")
            end
            
            hold off
            
        end
        
        
        function PlotState(obj, fig)
            
            if nargin > 1, figure(fig); hold on
            else, figure;
            end
            
            X = obj.State.XYZ.DR(1:2);
            cov = obj.State.Cov(1:2,1:2);
            
            plot( X(1), X(2), '*c');
            plotErrorEllipse(X, cov, 0.99, 'm');
            
        end
        
    end
end







function [theta, beta] = MinAngle(theta, beta)

if theta > beta, theta = theta - 2*pi;
else, beta = beta - 2*pi;
end

end





