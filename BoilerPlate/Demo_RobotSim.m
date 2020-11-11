 % Robot Sim
 % Jacob Anderson
 % 2/3/2020
 
close all
clear all
clc
format long g

ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(false);                                                     % Indicate whether new data should be selescted by the user


%% Get Data
bathFile = ui.GetFile('bathymetry','mat');                                  % Get file path GUI
if isempty(bathFile), return, end                                           % End script if there isn't an input for it to use

% ui = ui.NewData(true);    
wpFile = ui.GetFile('waypoints1','txt');                                   % Get file path GUI
if isempty(wpFile), return, end  

bathy = LoadBathymetry(bathFile);                                           % Load Bathymetry Map
bathy = bathy.ReduceResolution(0.5);


%% Create Robots

% ___ Vehicle kinematics __________________________________________________
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

F  = @(x,u,t,m) x + X(x,u,t) + M(m);                                        % Process Model

rad = @(x) x*pi/180;

% ___ Instanciate the AUV _________________________________________________
auv = R0b0t('AUV', 'auv1');                                                 % Instanciate the robot

auv = auv.SetDynamics('f', @(x,u,t,m) F(x, u, t, m) );                      % Add Process Model
auv = auv.SetDynamics('fx_dx', @(x,u,t,m) dX(x, u, t) );                    % Add Process Model derivative

auv = auv.SetUncertainty('Speed',   'Normal', -0.05, 0.5);                  % Set uncertainty in state transitions
auv = auv.SetUncertainty('Heading', 'Normal', rad(5), rad(4));

auv = auv.Add_Sensor("Depth",   'Normal', 0, 0.9);                          % Add Alitimeter with noise
auv = auv.Add_Sensor("Compass", 'Normal', rad(3), rad(6));                  % Add Compass with noise

auv.Nav = auv.Nav.Load_Waypoints(wpFile);                                   % Give the ASV a set of waypoints to follow

start = auv.Nav.Wp.Waypoints(1,1:2);                                        % Put the Robot at the First waypoint
auv   = auv.Set_VehicleState(start);

auv.Goal = auv.Nav.Wp.Waypoints(end,1:2);                                   % Put the Last waypoint as the goal  

clear bathFile wpFile d X dX M F H rad start                                % Clean up workspace


%% Diplay Initial Setup
fig = bathy.Plot_3DModel(0, 90);                                            % Show the Bathymetry Map
auv.Nav.PlotWaypoints(fig);
% auv.PlotLocation(fig);


%% Run Simulation
speed = 3;             % ASV1 Speed                        [m/s]

sim_time = 0;           % Simulation time                   [s]
dt = 0.5;               % Simulation time step              [s]
tic
while toc < 1000
    
    % ___ ASV Operations __________________________________________________
    auv = auv.Navigate(speed, dt , 3);
    
    % ____ Sensor Measurments _________________________________________
    auv = auv.Add_SensorMeasurment("Depth", -bathy.Depth( auv.State.XYZ.GT(1:2)'), sim_time);
    auv = auv.Add_SensorMeasurment("Compass",  auv.State.XYZ.GT(3), sim_time);
    
    % ____ Check for goal _________________________________________________
    if any( [auv.Goal_Achived(0.5), auv.Nav.Done]), break, end
    
    sim_time = sim_time +dt;                                                % time keeping
    
    % ____ Visulization ___________________________________________________
    if mod(sim_time, 50) == 0
        auv.PlotPaths(fig);
%         auv.Nav.PF.PlotState(fig, 0.5);
    end
    
%     auv.PlotLocation(fig);
     
end




%% Plotting
auv.PlotPaths(fig);



