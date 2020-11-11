% Demo ExtendedKalmanFilter
% Jacob Anderson
% 1/22/2020


clear all
close all
clc


%% System
d = 1.5;

X = @(x,u,t) [t*x(4)*cos(x(3));                                             % Process Model
              t*x(4)*sin(x(3)); 
              x(4)* tan(u(2))/d; 
                         t*u(1)];

dX = @(x,u,t) [ 1  0  -t*x(4)*sin(x(3))       t*cos(x(3));                  % Proces Model Jacobian
                0  1   t*x(4)*cos(x(3))       t*sin(x(3));
                0  0                  1  x(4)*tan(u(2))/d;
                0  0                  0                 1];    

M = @(m) [ m(2)*cos(m(1));                                                  % Process Noise Model
           m(2)*sin(m(1)); 
                        0; 
                        0];

dM = @(m) [-m(2)*sin(m(1))  cos(m(1));
            m(2)*cos(m(1))  sin(m(1));
                         0          0;
                         0          0];
                    
F  = @(x,u,t,m) x + X(x,u,t) + M(m);                                        % Process Model

H  = @(x,n) x + N(x) * n;                                                   % Sensor Model


%___ EKF __________________________________________________________________________________________________________________
ekf = ExtendedKalmanFilter('Demo');                                         % Instanciate the extended Kalman Filter

ekf = ekf.SetDynamics('f', @(x,u,t,m) F(x,u,t,m) );                         % Add Process Model
ekf = ekf.SetDynamics('h', @(x,n)     H(x,n) );                             % Add Sensor Model

ekf = ekf.SetNoise('f', 'Normal', 1, 4);

ekf = ekf.SetDerivative('df_dx', @(x,u,t) dX);                              % Derivaltive of process model with respect to x
ekf = ekf.SetDerivative('df_dm', @(x,u,t) dM(u));                            % Derivaltive of process model with respect to nois model
ekf = ekf.SetDerivative('dh_dx', @(x,u,t) 1 + N1(x));                       % Derivaltive of sensor model with respect to x
ekf = ekf.SetDerivative('dh_dn', @(x,u,t) N(x));                            % Derivaltive of sensor model with respect to noise model
%____________________________________________________________________________________________________________________________

clear M N N1

x   = [1;1;0;0];                                                                % State Vector
sig = [1, 1.5; 1.5, 3];                                                     % State Covariance
u   = [-2;-3];                                                              % Control Input
z   = [5; 5];                                                               % Observation

[x1, sig1] = ekf.Update(x, sig, u, z);                                        % Standared EKF Update

[x_b, sig_b] = ekf.BeliefUpdate(x, sig, u);                                     % EKF Belief Update






