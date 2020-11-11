% Autonomous Underwater Vehicle
% Jacob Anderson
% 2/11/2020


classdef AutonomousUnderwaterVehicle
    
    properties
        Name
        State
        Paths
        Sensors
        EKF
        Map
        VehicleProperties
        Nav
    end
    
    methods
%% Set Up
        
        % Constructor
        function obj = AutonomousUnderwaterVehicle(name)
            
            obj.Name = name;
            
            obj.VehicleProperties.fx    = @(x,u) x + u;                     % Default vehicle dynamics
            obj.VehicleProperties.fx_dx = @(x,u) 1;                         % Default vehicle dynamics
            
            obj.VehicleProperties.Limits.maxA = 2;                          % Max Acceleration
            obj.VehicleProperties.Limits.maxV = 5;                          % Max Velocity
            obj.VehicleProperties.Limits.maxTurn = 30 * pi/180;             % Max Turn Angle
            
            obj.VehicleProperties.Noise.Speed   = makedist('Normal','mu',0,'sigma',0.25);     % Default vehilce transition model
            obj.VehicleProperties.Noise.Heading = makedist('Normal','mu',0,'sigma',5);        % Default vehilce transition model
            %           obj.VehicleProperties.Noise.Pitch   = makedist('Normal','mu',0,'sigma',3);
            
            obj.State.XYZ.DR  = [0, 0, 0, 0]';                              % Deadreckoning Vehicle State [ x, y, heading, speed]
            obj.State.XYZ.GT  = [0, 0, 0, 0]';                              % Ground Trouth Vehicle State [ x, y, heading, speed]
            obj.State.Cov = zeros(4);
            
            obj.Paths.GT  = zeros(1000,4);
            obj.Paths.DR  = zeros(1000,4);
            obj.Paths.PF  = zeros(1000,4);
            obj.Paths.PF_indx  = 1;
            obj.Paths.ind = 1;
            
        end
        
        
        % Add Global Coordinates
        function obj = Add_Map(obj, map), obj.Map = map; end
        
        
        % Place Robt at certain place
        function obj = Set_VehicleState(obj, X)
            
            if numel(X) < 4, X = [X, 0, 0]; end
            
            obj.State.XYZ.DR = X';                                          % Set State
            obj.State.XYZ.GT = X';
            
            obj.Paths.DR(obj.Paths.ind,:) = obj.State.XYZ.DR;               % Show change in the path
            obj.Paths.GT(obj.Paths.ind,:) = obj.State.XYZ.GT;
        end
        
        
        % Span Robot in random location
        function obj = Spawn(obj, xlim, ylim)
            
            x = rand(4) .* [(max(xlim) - min(xlim)) + min(xlim);
                (max(ylim) - min(ylim)) + min(ylim);
                2*pi;
                obj.VehicleProperties.Limits.maxV];
            
            obj.State.XYZ.DR = x;
            obj.State.XYZ.GT = x;
            
        end
        
        
        % Function for use to specify system functions
        function obj = SetDynamics(obj, prop, func), obj.VehicleProperties.(prop) = func; end
        
        
        % Set vehicle transition model
        function obj = SetUncertainty(obj, type, dist, mu, sig)
            obj.VehicleProperties.Noise.(type) = makedist(dist,'mu',mu,'sigma',sig);
        end
        
        
        
%% Kinematicks
        
        % Move the vehicel towards a goal
        function obj = Move(obj, u, dt)
 
            % What the robot thinks is happening
            if isa(obj.EKF, 'ExtendedKalmanFilter')
                
                if isfield( obj.Sensors, 'GPS')
                    ii = obj.Sensors.GPS.iter;
                    z = obj.Sensors.GPS.value(ii,:)';
                    R = eye(4) * obj.Sensors.GPS.noise.sigma;
                    
                elseif isa( obj.Nav.PF, 'ParticleFilter')
                    z = obj.Nav.PF.X;
                    R = obj.Nav.PF.Cov;
                    
                else
                    z = zeros(size(obj.State.XYZ.DR));
                    R = Inf(numel(obj.State.XYZ.DR));
                end
                
                [X2, sig2] = obj.EKF.Update(obj.State.XYZ.DR, obj.State.Cov, u, z, dt, R);
                obj.State.XYZ.DR = X2;
                obj.State.Cov = sig2;
                
            else
                obj.State.XYZ.DR = obj.VehicleProperties.fx(obj.State.XYZ.DR, u, dt, [0,0]);
            end
            
            % What's Really Happening
            n = [ random(obj.Noise.Heading, 1), random(obj.Noise.Speed,1)];
            
            obj.State.XYZ.GT = obj.VehicleProperties.fx(obj.State.XYZ.GT, u, dt, n);
            
            
            
            % ____ Add Movement to path ___________________________________
            obj.Paths.DR(obj.Paths.ind,:) = obj.State.XYZ.DR;
            obj.Paths.GT(obj.Paths.ind,:) = obj.State.XYZ.GT;
            obj.Paths.ind = obj.Paths.ind + 1;
            
        end
        

        
%% Navigation Interface
        
        function obj = Navigate(obj, dt)
            
            
            if strcmpi(obj.Nav.Type, 'wp')
                
                x = obj.State.XYZ.DR(1:2)';
                
                [speed, bearing, obj.Nav.WP.indx, obj.Nav.Done] = WaypointFollow(x, obj.Nav.WP.Waypoints, obj.Nav.WP.Radius, obj.Nav.WP.indx);
                               
                % Signed shortes angle between the bearing and heading
                phi = (bearing - obj.State.XYZ.DR(3));
                
                if abs(phi) > 180
                    phi = mod( mod((phi),2*pi)+540 * pi/180, pi ) - pi;
                end
                
                if abs(phi) > obj.Limits.maxTurn
                    phi = obj.Limits.maxTurn * sign(phi);
                end
                
                % Control inputs
                a_ = (speed - obj.State.XYZ.DR(4)) / dt;
                a = min(abs(a_), obj.Limits.maxA);
                a = a * sign(a_);
                u = [a; phi];
                
            end
            
            obj = obj.Move(u, dt);
        end
        
        
        
        function obj = Load_Waypoints(obj, file, r)
            
            obj.Nav.Type = 'wp';
            obj.Nav.WP.Waypoints = csvread(file);
            obj.Nav.WP.Radius = r;
            obj.Nav.WP.indx = 1;
            obj.Nav.Done = false;
            
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
        

        
%% Plotting 
               
% Plot the location of the robot
function PlotLocation(obj, fig)
    
    figure(fig)
    hold on
    
    
    s1 = scatter3(obj.State.XYZ.DR(1), obj.State.XYZ.DR(2), obj.State.XYZ.DR(3), '.c', ...
        'MarkerEdgeAlpha', 0.6, ...
        'MarkerFaceAlpha', 0.6);
    s2 = scatter3(obj.State.XYZ.GT(1), obj.State.XYZ.GT(2), obj.State.XYZ.GT(3), '.r');
    legend([s1,s2], 'DR', 'GT')
    
    
    
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
                    p2 = plot3(obj.Paths.DR(1:i,1), obj.Paths.DR(1:i,2), ones(i,1), 'r'); %obj.Paths.DR(1:i,3), 'r')
                    p3 = plot3(obj.Paths.PF(1:i,1), obj.Paths.PF(1:i,2), ones(i,1), 'b'); %obj.Paths.DR(1:i,3), 'b')
                    legend([p1, p2, p3], 'GT', 'DR', 'PF')
                    
                otherwise
                    warning("We Can Not Disclose the locaiton of this device")
            end
            
            hold off
            
        end
        
        
            
    end
end