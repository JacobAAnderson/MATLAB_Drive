 % Robot with Particle Filter and EKF localization Sim
 % Jacob Anderson
 % 1/30/2020
 
close all
clear all
clc
format long g

ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(true);                                                     % Indicate whether new data should be selescted by the user


%% Get Data
bathFile = ui.GetFile('bathymetry','mat');                                  % Get file path GUI
if isempty(bathFile), return, end                                           % End script if there isn't an input for it to use

% ui = ui.NewData(true);    
wpFile = ui.GetFile('waypoints1','txt');                                    % Get file path GUI
if isempty(wpFile), return, end  

bathy = LoadBathymetry(bathFile);                                           % Load Bathymetry Map
bathy = bathy.ReduceResolution(0.5);

ti = Terrain_Info(bathy.AsMapp);
ti = ti.Eval_MapInfo(20);

%% Create Robot

auv = AUV1( 'AUV', bathy, wpFile);

rad = @(x) x*pi/180;

%% Diplay Initial Setup
fig = bathy.Plot_3DModel(0, 90);                                            % Show the Bathymetry Map
auv.Nav.PlotWaypoints(fig);
% auv.PlotLocation(fig);
hold on

sig_s = 0.3;
sig_h = rad(5);
path = auv.Nav.Wp.Waypoints(:,1:2);
wp_idx = 1;
X = path(1,:);
Cov = zeros(2);
sig_alt = 0.5;

logistics = @(x) 2/(1+exp(-6*x+1))-1;


%% Run Simulation
auv_speed = 3;             % ASV1 Speed                        [m/s]

sim_time = 0;           % Simulation time                   [s]
dt = 0.5;               % Simulation time step              [s]
tic
while toc < 1000
    
    % ____ Sensor Measurments _________________________________________
    auv = auv.Add_SensorMeasurment("Depth",   -bathy.Depth( auv.State.XYZ.GT(1:2)'), sim_time);
    auv = auv.Add_SensorMeasurment("Compass", auv.State.XYZ.GT(3), sim_time);
    
    
    % ___ Update Particle Filter __________________________________________
    auv = auv.UpdatePF(dt);
    
    % ___ ASV Operations __________________________________________________
    auv = auv.Navigate(auv_speed, dt , 3);
    
    
    % ____ Check for goal _________________________________________________
    if any( [auv.Goal_Achived(0.5), auv.Nav.Done]), break, end
    
    sim_time = sim_time +dt;                                                % time keeping
    
    % ____ Visulization ___________________________________________________
    if mod(sim_time, 50) == 0
        auv.PlotPaths(fig);
        auv.Nav.PF.PlotState(fig, 0.5);
        
        hold on
        plotErrorEllipse(X, Cov, 0.99, 'y')
        drawnow
    end
    
%     auv.PlotLocation(fig);
     

    sig_bathy = ti.Map.Measure(X, "STD");

    [speed, heading, wp_idx, ~] = WaypointFollow(X', path, 3, wp_idx);
    X = DeadReckoning(X, speed * 3, heading, dt);
    Cov = Cov + logistics(sig_alt/sig_bathy) * Est_Cov(speed, sig_s, heading, sig_h, dt);
    
    plot(X(1),X(2),'.m')
end




%% Plotting
auv.PlotPaths(fig);






