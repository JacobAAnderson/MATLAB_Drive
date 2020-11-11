classdef ExtendedKalmanFilter
    
    properties
        name
%     end
%     
%     properties (Access = private)
        dynamics
        noise
        derivatives
    end
    
    methods
        
        function obj = ExtendedKalmanFilter(name)
            
            obj.name = name;
            
            obj.noise.f = makedist('Normal','mu',0,'sigma',1);              % Instanciate default values
            obj.noise.h = makedist('Normal','mu',0,'sigma',1);
            
            obj.dynamics.f = @(x,u) x + u + random(obj.noise.f,1,2);
            obj.dynamics.h = @(x) x + random(obj.noise.h,1,2);
            
            obj.derivatives.df_dx = @(x,u) 1;
            obj.derivatives.df_dm = @(x,u) 1;
            
            obj.derivatives.dh_dx = @(x,u) 1;
            obj.derivatives.dh_dn = @(x,u) 1;
            
        end
        
        
        function obj = SetDynamics(obj, prop, func), obj.dynamics.(prop) = func; end    % Function for use to specify system functions
        
        
        function obj = SetDerivative(obj, prop, func), obj.derivatives.(prop) = func; end   % Set Derivatives for determining kalman gain
        
        
        function obj = SetNoise(obj, prop, dist, mu, sig)                   % Set Noise Distributions
            obj.noise.(prop) = makedist(dist,'mu',mu,'sigma',sig);
        end
        
        
        function [x_, sig_] = Update(obj, x, sigma, u, z, dt, R)
            
            
%             fprintf('\n\nR:')
%             disp(R)
%             fprintf('\n\n')
            
            if nargin < 6, dt = 1; end
            
            m = [0;0];
            
            A = obj.derivatives.df_dx(x,u,dt);                            % Evaluate Derivative terms
            M = obj.derivatives.df_dm(x,u,dt);
            H = obj.derivatives.dh_dx(x,u,dt);
            N = obj.derivatives.dh_dn(x,u,dt);
            
            P = A * sigma * A' + M * M';                                    % Covariance Term
            
            if nargin < 7
                K = P * H' * pinv(H * P * H' + N * N');                     % Kalman Gain
            else
                K = P * H' * pinv(H * P * H' + R);                          % Kalman Gain
            end
            
            x_ = obj.dynamics.f(x, u, dt, m);                               % Updated State
            x_ = x_ + K*(z - obj.dynamics.h(x_, x_ .*0));
            
            sig_ = (eye(size(K)) - K*H) *P;                                 % New Covariance
            sig_ = (sig_ + sig_')/2;                                        % Make sure the matrix is symmetric positive semi-definite
            
        end
        
        
        function [x_, sig_] = BeliefUpdate(obj, x, sigma, u, dt)
            
            if nargin < 5, dt = 1; end
            
            A = obj.derivatives.df_dx(x,u,dt);                              % Evaluate Derivative terms
            M = obj.derivatives.df_dm(x,u,dt);
            H = obj.derivatives.dh_dx(x,u,dt);
            N = obj.derivatives.dh_dn(x,u,dt);
            
            P = A * sigma * A' + M * M';                                    % Covariance Term
            K = P * H' * pinv(H * P * H' + N * N');                         % Kalman Gain
            
            W = K*H*P;                                                      % Noise Term
            W = ( W + W')/2;                                                % Make sure the matrix is symmetric positive semi-definite
            R = mvnrnd([0,0],W,1);                                          % Generate Noise
            
            x_ = obj.dynamics.f(x, u, dt, 0) + R';                          % Updated State
            
            sig_ = (eye(size(K)) - K*H) *P;                                 % New Covariance
            sig_ = (sig_ + sig_')/2;                                        % Make sure the matrix is symmetric positive semi-definite
            
        end
        
    end
    
end



