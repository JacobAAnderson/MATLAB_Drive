% Demo ExtendedKalmanFilter
% Jacob Anderson
% 1/22/2020


clear all
close all
clc


M  = @(u) [0.5, 0.1; 0.1, 0.5] * u;                                         % Process Noise Model

N  = @(x) [1./(1+exp(5-x(1))), 0;                                           % Sensor Noise  Model
           0, 1./(1+exp(5-x(1)))];

N1 = @(x) [exp(5 - x(1)) ./ ((1+exp(5-x(1))) .* (1+exp(5-x(1))) ), 0;       % 1st derivative of sensor Noise Model
           0, exp(5 - x(1)) ./ ((1+exp(5-x(1))) .* (1+exp(5-x(1))) )];

F  = @(x,u,t,m) x + u + M(u).*m;                                            % Process Model
H  = @(x,n) x + N(x) * n;                                                   % Sensor Model


%___ EKF __________________________________________________________________________________________________________________
ekf = ExtendedKalmanFilter('Demo');                                         % Instanciate the extended Kalman Filter

ekf = ekf.SetDynamics('f', @(x,u,t,m) F(x,u,t,m) );                         % Add Process Model
ekf = ekf.SetDynamics('h', @(x,n)     H(x,n) );                             % Add Sensor Model

ekf = ekf.SetNoise('f', 'Normal', 1, 4);

ekf = ekf.SetDerivative('df_dx', @(x,u,t) eye(2));                          % Derivaltive of process model with respect to x
ekf = ekf.SetDerivative('df_dm', @(x,u,t) M(u));                            % Derivaltive of process model with respect to nois model
ekf = ekf.SetDerivative('dh_dx', @(x,u,t) 1 + N1(x));                       % Derivaltive of sensor model with respect to x
ekf = ekf.SetDerivative('dh_dn', @(x,u,t) N(x));                            % Derivaltive of sensor model with respect to noise model
%____________________________________________________________________________________________________________________________

clear M N N1

x   = [1;1];                                                                % State Vector
sig = [1, 1.5; 1.5, 3];                                                     % State Covariance
u   = [-2;-3];                                                              % Control Input
z   = [5; 5];                                                               % Observation

[x1, sig1] = ekf.Update(x, sig, u, z);                                        % Standared EKF Update

[x_b, sig_b] = ekf.BeliefUpdate(x, sig, u);                                     % EKF Belief Update






