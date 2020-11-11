classdef ParticleFilter00
    
    properties
        stats                                                               % Statistics on the particles
        Cov
    end
    
    properties (Access = private)
        particles                                                           % List of partilce [easting (utm), northing (utm), weight]
        map                                                                 % Structure containing the map
        num_Particles                                                       % Number of particles
        noise                                                               % Vehicle transition Model as a Gaussian
        iter
    end
    
    
    methods
        
        % Initilaize the Partilce filter
        function obj = ParticleFilter00( num_Particles, map, iters )
            
            obj.map = map;                                                  % Map structure
            obj.num_Particles = num_Particles;                              % Number of Partilces
            
            % Randomly Initilaze the particles location accross the map
            obj.particles(:,1:3) = rand(obj.num_Particles, 3) .* [(obj.map.uperRight - obj.map.origin), 0] + [obj.map.origin, 0];
            obj.particles(:,3) = 1/num_Particles;                           % Normalixe the weights' distribuiton
            
            % Models
            obj.noise.sensor  = makedist('Normal','mu',0,'sigma',0.25);     % Sensor Model default values
            obj.noise.vehicle = makedist('Normal','mu',0,'sigma',0.25);     % Default vehilce transition model
            obj.noise.compass = makedist('Normal','mu',0,'sigma',3);        % Default vehilce transition model
            
            if nargin > 2
                obj.iter = 1;
                obj.stats.mean    = zeros(iters, 3);
                obj.stats.std     = zeros(iters, 3);
                obj.stats.max_min = zeros(iters, 2);
                obj.Cov = cov(obj.particles(:,1), obj.particles(:,2));
            else
                obj.iter = NaN;
                obj.stats = NaN;
            end

        end
        
        
        % Initialize the particles' location
        function obj = SetLocation(obj, location)
            obj.particles(:,1) = location(1);
            obj.particles(:,2) = location(2);
            obj.particles(:,3) = 1/obj.num_Particles;
        end
        
        
        % Adjust State estimate
        function obj = AdjustLocation(obj, location)
           
            d_xy = location.pose - obj.GetLocation;
            
            obj.particles(:, 1:2) = obj.particles(:, 1:2) + d_xy;
            obj.particles(:,3) = 1/obj.num_Particles;                           % Normalixe the weights' distribuiton
            
        end
               
        
        % Set the values of the sensor model
        function obj = SetNoise(obj, type, dist, mu, sig)
            obj.noise.(type) = makedist(dist,'mu',mu,'sigma',sig);
        end
  
        
        % Update Partilce Filter Locations and Weigths
        function obj = Update(obj, velocity, heading, dT, measurment, offSet)
            
            % Move the particles
            heading = heading + random(obj.noise.compass, obj.num_Particles, 1);
            heading(heading > 360) = heading(heading > 360) - 360;
            heading(heading < 0)   = heading(heading < 0)   + 360;
            
            dXdY = velocity * dT * [sind(heading), cosd(heading)];          % Dead reckoning movement
            
            noise_ = random(obj.noise.vehicle, obj.num_Particles, 2);       % Add noise based on the vehicle's transition model
            
            newParticles(:,1:2) = obj.particles(:, 1:2) + dXdY + noise_;
            
            % Prevent Particles from leaving the map
            newParticles(newParticles(:,1) > obj.map.uperRight(1), 1) = obj.map.uperRight(1);
            newParticles(newParticles(:,2) > obj.map.uperRight(2), 2) = obj.map.uperRight(2);
            
            newParticles(newParticles(:,1) < obj.map.origin(1), 1) = obj.map.origin(1);
            newParticles(newParticles(:,2) < obj.map.origin(2), 2) = obj.map.origin(2);

            % Sample the map
            samples = newParticles + offSet;                                % Account for vehicle attitude
            depths  = obj.map.Depth(samples);
            
            % Update the weights
            diffs = depths - abs(measurment);
            
            if any(isnan(diffs))
                disp('NaNs')
            end
            
            weights = pdf(obj.noise.sensor, diffs);                         % Assign weights based on sensor model
            weights(isnan(weights)) = 0;                                    % Get Rid of nan values
            
            if all( weights == 0)                                           % Bad sensor measurement creates 0 probabilities
                obj.particles = [newParticles, obj.particles(:, 3)];        % Retain Previus Weights
                
            else                                                            % Good-ish sensor measuremnt creates normal probailities
                weights = weights / sum(weights);                           % Normalize the weigths
                obj.particles = [newParticles, obj.particles(:, 3) .* weights];
            end
            
            
            if (max(heading) - min(heading)) > 180
                heading(heading< 180) = heading(heading<180) + 360;
            end
            
            heading = heading * 2*pi / 360;
            
            obj.Cov = cov([obj.particles(:,1), obj.particles(:,2),heading]); % Calculate the covariance of the particles
            
            if ~isnan(obj.iter)                                             % Keep track of stats
                obj.iter = obj.iter+1;
                obj.stats.mean(obj.iter,:) = [nanmean(obj.particles(:,1)), ...
                                              nanmean(obj.particles(:,2)), ...
                                              nanmean(obj.particles(:,3))];
                                          
                obj.stats.std(obj.iter,:) = [nanstd(obj.particles(:,1)), ...  .*obj.particles(:,3)), ...  
                                             nanstd(obj.particles(:,2)), ...  .*obj.particles(:,3)), ...  
                                             nanstd(obj.particles(:,3))];
                                         
                obj.stats.max_min(obj.iter,:) = [nanmax( obj.particles(:,3) ), nanmin( obj.particles(:,3) )];
            end
            
            
            obj = obj.ReSample;                                             % Resample to get rid of 0s
            
            
        end
        
        
        % Get Vehicle Location
        function location = GetLocation(obj)
            
            location = [nansum(obj.particles(:,1) .* obj.particles(:,3)), nansum(obj.particles(:,2) .* obj.particles(:,3))];
            
        end
        
        
        % Plot Particle filter
        function fig = Plot(obj, emData, plotView)
            
            location = obj.GetLocation;                                     % Get vehicle location
            
%             isgraphics(obj.fig,'figure')
            
            % Plot Bathymetry Map
            if nargin < 3
                obj.fig = obj.map.Plot_3DModel(250, 50);
            else
                obj.fig = obj.map.Plot_3DModel(plotView(1), plotView(2));
            end
            
            
            hold on
            plot( obj.particles(:,1), obj.particles(:,2), '.r')             % Plot Particles
            plot( location(1), location(2), '+b')                           % Plot Particle filter's estimate of the vehicle's location
            plot( location(1), location(2), 'ob')
            
            % Plot Ground truth data if available
            if nargin > 1
                plot(emData(:,1), emData(:,2))
                plot(emData(end,1), emData(end,2),'+g')
                plot(emData(end,1), emData(end,2),'og')
            end
            
            hold off
            drawnow
            
            if nargout > 0
                fig = obj.fig;
            end
            
        end

        
        function fig = PlotState(obj, fig, trans)
            
            if nargin > 1, figure(fig), hold on;
            else, fig = figure('name',"Particle Filter State");
            end
            
            s = scatter(obj.particles(:,1), obj.particles(:,2), '.w');
            
            if nargin > 2
                s.MarkerEdgeAlpha = trans;
                s.MarkerFaceAlpha = trans;
            end
            
            hold on
            
            xy = obj.GetLocation;
            plot(xy(1), xy(2), '+r');
            plotErrorEllipse(obj.GetLocation, obj.Cov(1:2,1:2), 50, 'm');
            drawnow
            
            hold off
        end
            
        
    end
    
    
    methods (Access = private)
        
        % Resample the particles
        function obj = ReSample(obj)
            
            threshold = mean(obj.particles(:,3)) * 0.001 / obj.num_Particles; % Threshold for which particles to keep
            
            in = obj.particles(:,3) >= threshold;                           % Keep particles with higher weights
            
            keep = obj.particles( in, :);                                   % Particles to keep
            
            if sum(in) < min(1, obj.num_Particles/100)
                
%                 location = obj.GetLocation;  
                
%                 newParticles(:,1) = rand(obj.num_Particles/2, 1) * (obj.uperRight(1)-obj.origin(1)) + obj.origin(1);
%                 newParticles(:,2) = rand(obj.num_Particles/2, 1) * (obj.uperRight(2)-obj.origin(2)) + obj.origin(2);
%                 newParticles(:,3) = 1/obj.num_Particles;                           % Normalixe the weights' distribuiton
                
                keep = [keep; [obj.GetLocation, 1/obj.num_Particles]];

            end
            
            
            
            n = ceil( obj.num_Particles / size(keep,1));                    % How many times the remaning particles need to be reproduced to maintina the size of the particle filter
            new = repmat(keep, n,1);                                        % Make new particles
            keep = [keep; new];
            keep = keep(1:obj.num_Particles, :);
            
            keep(:,3) = keep(:,3) / nansum(keep(:,3));                      % Normalize the weights
            
            obj.particles = keep;
           
        end
        
    end
    
    
    
end