% Belief Space POMDP
% Jacob Anderson
% 01/14/2020

close all
clear all
clc


%% Create Environment
x = 0:0.01:20;
y = 0:0.01:10;

[X,Y] = meshgrid(x,y);

Z = 1./(1+exp(5-X));


wall1 = [14, 10; 14, 5.75; 15, 5.75; 15,10];
wall2 = [14,  0; 14, 4.25; 15, 4.25; 15, 0];

map = Mapp(X, Y, Z, 'Environment');
map = map.Add_Obsticle(wall1);
map = map.Add_Obsticle(wall2);

map.PlotMap;

clear x X y Y Z wall1 wall2

[hit, dev] = map.Collition([10,5], [9,6;6,9]);



%% Nominal Path

start = [11.2673, 9.0816];
goal  = [17.8571, 1.4723];

rrt = RRT(map, 3);

path = rrt.Plan(start, goal, 1000, false);

hold on
plot(start(1), start(2), '*m')
plot(goal(1), goal(2), '*g')
plot(path(:,1), path(:,2), 'b')
hold off

% OR --------------------------------------------

path = [11.2673, 9.0816;
        12.6037, 8.7318;
        12.3733, 7.7697;
        11.5438, 6.9825;
        12.5576, 6.1953;
        13.9401, 5.0292;
        15.7373, 4.6793;
        16.7972, 3.8630;
        17.0276, 2.4344;
        17.8571, 1.4723];

hold on
plot(path(1,1), path(1,2), '*r')
plot(path(end,1), path(end,2), '*g')
plot(path(:,1), path(:,2))
hold off




%% System Functions
M  = @(u) [0.5, 0.1; 0.1, 0.5] * u;                                         % Process Noise Model
N  = @(x) [1./(1+exp(5-x(1))), 0;                                           % Sensor Noise  Model
           0, 1./(1+exp(5-x(1)))];
N1 = @(x) [exp(5 - x(1)) ./ ((1+exp(5-x(1))) .* (1+exp(5-x(1))) ), 0;       % 1st derivative of sensor Noise Model
           0, exp(5 - x(1)) ./ ((1+exp(5-x(1))) .* (1+exp(5-x(1))) )];

F  = @(x,u, m) x + u + M(u)*m;                                              % Process Model
H  = @(x, n) x + N(x) * n;                                                  % Sensor Model


%___ EKF __________________________________________________________________________________________________________________
ekf = ExtendedKalmanFilter('Belief');                                       % Instanciate the extended Kalman Filter

ekf = ekf.SetDynamics('f', @(x,u,m) F(x, u, m) );                           % Add Process Model
ekf = ekf.SetDynamics('h', @(x, n) H(x, n) );                               % Add Sensor Model

ekf = ekf.SetDerivative('df_dx', @(x,u) eye(2));                            % Derivaltive of process model with respect to x
ekf = ekf.SetDerivative('df_dm', @(x,u) M(u));                              % Derivaltive of process model with respect to nois model
ekf = ekf.SetDerivative('dh_dx', @(x,u) 1 + N1(x));                         % Derivaltive of sensor model with respect to x
ekf = ekf.SetDerivative('dh_dn', @(x,u) N(x));                              % Derivaltive of sensor model with respect to noise model
%____________________________________________________________________________________________________________________________

clear M N N1

x   = [1;1];                                                                % State Vector
sig = [1, 1.5; 1.5, 3];                                                     % State Covariance
u   = [-2;-3];                                                              % Control Input
z   = [5; 5];                                                               % Observation

[x, sig] = ekf.Update(x, sig, u, z);                                        % Standared EKF Update

[x, sig] = ekf.BeliefUpdate(x, sig, u);                                     % EKF Belief Update








%% Value Iteration




%% Functions

function c = C_l(x,sig)
q = [0.5,0.5];
c = x'* q * x + trace(q*sig);
end
    
    
function c = C_t(x, sig, u, map)
r = [0.5,0.5];
q = [0.5,0.5];

[~, dev] = map.Collition(x, sig);

Y = gammainc(numel(x)/2, dev*dev/2,'lower') / gamma(numel(x)/2); 

h = -log10(Y);

c = u'*r*u + trace(q*sig) + h;

end

