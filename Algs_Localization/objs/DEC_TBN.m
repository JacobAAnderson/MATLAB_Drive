classdef DEC_TBN
    
    properties
        cov
        X
        name
        particles
        Map
        path
        path_idx
        noise
    end
    
    properties (Access = private)
        num_Particles
        
    end
    
    
    methods
        
        % Initilaize the Partilce filter
        function obj = DEC_TBN(name, num_Particles)
            
            obj.name = name;
            obj.num_Particles = num_Particles;                              % Number of Partilces
            
            % Models
            obj.noise.speed     = makedist('Normal','mu',0,'sigma',0.25);       % Default speed model
            obj.noise.compass   = makedist('Normal','mu',0,'sigma',3);          % Default compass model
            obj.noise.altimeter = makedist('Normal', 'mu', 0, 'sigma', 0.25);   % Default altimiter model
            
            % Create Particles
            obj.particles(:,1:3) = zeros(obj.num_Particles, 3);
            obj.particles(:,3) = 1/num_Particles;                           % Normalixe the weights' distribuiton
            
            % Create Path
            obj.path = zeros(500,4);
            obj.path_idx = 1;
            
        end
        
        
        function obj = Add_Map(obj,map)                                     % Add Map and initialize the particle filter
            obj.Map = map;
            obj.particles(:,1:2) = rand(obj.num_Particles, 2) .* (map.uperRight - map.origin) + map.origin;
        end
        
         
        function obj = SetNoise(obj, prop, dist, mu, sig)                   % Set Custom Noise Distributions
            obj.noise.(prop) = makedist(dist,'mu',mu,'sigma',sig);
        end
        
        
        % Initialize the particles' location
        function obj = SetLocation(obj, location, sig)
            
            obj.particles = [mvnrnd(location,sig,obj.num_Particles,1), ones(obj.num_Particles,1) /obj.num_Particles];
            obj.X = [location(1);location(2); 0; 0];
            obj.cov = [sig, zeros(2); zeros(2), eye(2) * 0.25];
            
        end
        
        
        % Adjust State estimate
        function obj = AdjustLocation(obj, location)
            
            d_xy = location.pose - obj.GetLocation;
            
            obj.particles(:, 1:2) = obj.particles(:, 1:2) + d_xy;
            
            obj.particles(:,3) = 1/obj.num_Particles;
            
        end
        
        
        % Update Partilce Filter Locations and Weigths from Dead Reckoning
        function obj = Update(obj, speed, dt, heading, measurment, offset)
            
            % --- Move the particles ---
            heading = heading + random(obj.noise.compass, obj.num_Particles, 1); % Add noise to the heading
            heading(heading > 360) = heading(heading > 360) - 360;               % Make Sure headings are between 0 and 360
            heading(heading < 0)   = heading(heading < 0)   + 360;
            
            speed = speed + random(obj.noise.speed, obj.num_Particles, 1);       % Add noise to the speed
            
            newParticles = DeadReckoning(obj.particles(:, 1:2), speed, heading, dt, 'deg', 'compass');
            
            % Prevent Particles from leaving the map
            newParticles(newParticles(:,1) > obj.Map.uperRight(1), 1) = obj.Map.uperRight(1);
            newParticles(newParticles(:,2) > obj.Map.uperRight(2), 2) = obj.Map.uperRight(2);
            
            newParticles(newParticles(:,1) < obj.Map.origin(1), 1) = obj.Map.origin(1);
            newParticles(newParticles(:,2) < obj.Map.origin(2), 2) = obj.Map.origin(2);
            
            % --- Update weights based on altimeter reading if provided ---
            if nargin >= 5
                
                if nargin == 6
                    depths  = obj.Map.Depth(newParticles + offset);         % Sample the map, apply attitude offset
                else
                    depths  = obj.Map.Depth(newParticles);                  % Sample the map
                end
                
                diffs = depths + measurment;                                % Update the weights
                
%                 if any(isnan(diffs))
%                     warning('NaN measurments from bathimetry map')
%                 end
                
                weights = pdf(obj.noise.altimeter, diffs);                  % Assign weights based on sensor model
                weights(isnan(weights)) = 0;                                % Get Rid of nan values
                
                if all( weights == 0)                                       % Bad sensor measurement creates 0 probabilities
                    obj.particles(:,1:2) = newParticles;                      % Retain Previus Weights
                    
                else                                                        % Good-ish sensor measuremnt creates normal probailities
                    weights = obj.particles(:, 3) .* weights;               % Update Particle weights
                    weights = weights / sum(weights);                       % Normalize the weigths
                    obj.particles = [newParticles, weights];
                end
            
                
            else % Otherwise, just move the particles
                obj.particles(:,1:2) = newParticles;                        
            end
            
            
            [obj.X, obj.cov] = Covariance(obj.particles, heading, speed);
            
            Neff = 1/sum(obj.particles(:, 3) .* obj.particles(:, 3));       % Calculate the number of effective particles
            Nt = 0.50 * obj.num_Particles;
            
            if Neff < Nt, obj = obj.ReSample; end                           % Resample if the number of effective particles is bellow the threshold
        
            obj.path(obj.path_idx,:) = obj.X;                               % Keep track of path
            obj.path_idx = obj.path_idx + 1;
            
        end
        
        
        
        % Update Partilce Filter Weigths from Acomms communication
        function obj = Acoms(obj, dist, x, cov_)
            
            particles_ = obj.particles(:, 1:2);
            weights    = obj.particles(:, 3);                                % Assign weights based on sensor model
            
            % Acoms Update
%             dist = acoms{2};
%             x    = acoms{3};
%             cov_ = acoms{4};
            
            distances = x - particles_;                                     % Project the range measurments onto the new particles
            distances = distances./ vecnorm(distances, 2, 2);
            samplePoints = particles_ + distances * dist;
            
            s = obj.noise.acoms.sigma;                                      % Include acoms error
            s = [s*distances(1), 0; 0, s*distances(2)];
           
            cov_ = sqrt( cov_*cov_ + s*s);                                  % Make sure the matrix is symmetric positive semi-definite
            cov_ = (cov_ + cov_')/2;      
            
            W = mvnpdf(samplePoints, x, cov_);                              % Likelyhood of the acoms measurments
            
            weights = weights .* W;
            weights(isnan(weights)) = 0;                                    % Get Rid of nan values
            
            if all( weights == 0), return, end                              % Bad sensor measurement creates 0 probabilities, do not update particles
 
            
            weights = obj.particles(:, 3) .* weights;                       % Update Particle weights
            weights = weights / sum(weights);                               % Normalize the weigths
            obj.particles = [particles_, weights];
            
            Neff = 1/sum(obj.particles(:, 3) .* obj.particles(:, 3));
            Nt = 0.50 * obj.num_Particles;
            
            if Neff < Nt, obj = obj.ReSample; end                           % Resample to get rid of 0s
            
            
            
        end
        
        
        % Plot Particle filter
        function fig = Plot(obj, emData, plotView)
            
            location = obj.X;                                               % Get vehicle location
            
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
            
            plot(obj.X(1), obj.X(2), '+r');
            plotErrorEllipse(obj.X, obj.cov(1:2,1:2), 0.99, 'm');
            drawnow
            
            hold off
        end
        
        
        function fig = PlotPath(obj, fig)
            
            if nargin > 1, figure(fig), hold on;
            else, fig = figure('name',"Particle Filter Path", 'numbertitle','off');
            end
            
            path_ = obj.path(1:obj.path_idx-1, 1:2);
            
            plot(path_(:,1), path_(:,2))
            
            hold off
            
            
        end
        
    end
    
    
    methods (Access = private)
        
        % Resample the particles
        function obj = ReSample(obj)
           
            threshold = mean(obj.particles(:,3)) * 0.001 / obj.num_Particles; % Threshold for which particles to keep
            
            in = obj.particles(:,3) >= threshold;                           % Keep particles with higher weights
           
%             out = sum(~in);            
%             obj.particles(~in,:) = [mvnrnd(obj.X(1:2), 3.* obj.cov(1:2,1:2), out,1), ones(out,1)/obj.num_Particles];
%             obj.particles(:,3) = 1/obj.num_Particles;
            
            keep = obj.particles( in, :);                                   % Particles to keep
                        
            if sum(in) < min(1, obj.num_Particles/100)
                
%                 location = obj.GetLocation;
%                 
%                 newParticles(:,1) = rand(obj.num_Particles/2, 1) * (obj.uperRight(1)-obj.origin(1)) + obj.origin(1);
%                 newParticles(:,2) = rand(obj.num_Particles/2, 1) * (obj.uperRight(2)-obj.origin(2)) + obj.origin(2);
%                 newParticles(:,3) = 1/obj.num_Particles;                           % Normalixe the weights' distribuiton
                
                keep = [keep; [obj.X(1:2)', 1/obj.num_Particles]];
                
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




function [X, Cov] = Covariance(particles, headings, speeds)

% Calculate the covariance of the particles
if (max(headings) - min(headings)) > 180                                    % Make sure the headings do not cross 360/0 when calculating covariance
    headings(headings< 180) = headings(headings<180) + 360;
end

theta = nansum(headings .* particles(:,3));                                 % Get heading with out 360/0 crossovers

% t = theta;                                                                  % Save this heading for covariance calculations

theta(theta > 360) = theta(theta > 360) - 360;                              % Make sure heading is between 0 and 360
theta(theta < 0)   = theta(theta < 0)   + 360;

X = [nansum(particles(:,1) .* particles(:,3));      % X
     nansum(particles(:,2) .* particles(:,3));      % Y
                                        theta;      % Heading as a compass angle in degrees
     nansum(        speeds .* particles(:,3))];     % Speed


particles = [particles(:,1:2), headings *pi/180, speeds]; 
Cov = cov(particles);

% x = X;
% x(3) = t *pi/180;

% w = particles(:,3);
% particles = [particles(:,1:2), headings *pi/180, speeds] - x';

% Cov = zeros(4);
% for ii = 1: size(particles,1)
%     Cov = Cov + w(ii) .* (particles(ii,:)' * particles(ii,:));
% end

end

