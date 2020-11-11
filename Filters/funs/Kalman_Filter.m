% Kalman Filter
% Jacob Anderson
% 2/7/2020




 function [X, Sigma] = Kalman_Filter(x, sig, u, z, A, B, C, Q, R)
           
           X_ = A * x  +  B * u;
           P_ = A * sig * A'  +  Q;
           
           % K = Err in estimate / ( err in eatimate + err in measurment)
           K = P_ * C' * pinv(C * P_ * C' + R);
           
           X = X_ + K * (z - C * X_);
           
           % sigma = (1-K) * Err in estimate
           Sigma = (eye(size(A)) - K * C) * P_;
           
           
       end



%% classdef Kalman_Filter
%     
%    properties
%        Name
%        A                % Process Model
%        B                % Control Matrix
%        C                % Observation Matrix
%        Q                % Process Standard Deviation
%        R                % Measurment Standard Deviation
%    end
%    
%    methods
%        
%        function obj = Kalman_Filter(a, b, c, name)
%            obj.A = a;
%            obj.B = b;
%            obj.C = c;
%            
%            if nargin > 3, obj.Name = name; end
%            
%        end
%     
%        function [X, Sigma] = Update(obj, x, sig, u, z, R)
%            
%            if nargin < 6, R = obj.R; end
%            
%            X_ = obj.A*x + obj.B*u;
%            P_ = obj.A*sig*obj.A' + obj.Q;
%            
%            % K = Err in estimate / ( err in eatimate + err in measurment)
%            K = P_ * obj.C' * pinv(obj.C * P_ * obj.C' + R);
%            
%            X = X_ + K * (z - obj.C * X_);
%            
%            % sigma = (1-K) * Err in estimate
%            Sigma = (eye(size(obj.A)) - K * obj.C) * P_;
%            
%            
%        end
%     
%     
%    end
% end