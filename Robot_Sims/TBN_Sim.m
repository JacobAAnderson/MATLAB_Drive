% Robot Simulator

close all
clear all
clc
format long g

ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(false);                                                      % Indicate whether new data should be selescted by the user

bathFile = ui.GetFile('bathymetry','mat');                                  % Get file path GUI
if isempty(bathFile), return, end                                           % End script if there isn't an input for it to use

% ui = ui.NewData(true);    
wpFile = ui.GetFile('waypoints','txt');                                     % Get file path GUI
if isempty(wpFile), return, end  

bathy = LoadBathymetry(bathFile);                                           % Load Bathymetry Map
clear bathFile                                                              % Clean up workspace

bathy = bathy.ReduceResolution(0.5);

% Make terrain info maps
ti = Terrain_Info(bathy.AsMapp);
ti = ti.Eval_MapInfo(20);


%% Create ASV
asv = R0b0t('ASV');

asv = asv.SetUncertainty('Speed',   'Normal', 0, 0.25);
asv = asv.SetUncertainty('Compass','Normal', 0, 10);

asv = asv.LoadWaypoints(wpFile);                                            % Give the ASV a set of waypoints to follow
start = asv.Nav.Wp.Waypoints(1,:);                                          % Put th eRobot as the First waypoint
 
asv = asv.SetLocation(start(1), start(2));

asv = asv.Add_Sensor("Depth",   'Normal', 0, 0);
asv = asv.Add_Sensor("Entropy", 'Normal', 0, 0);
asv = asv.Add_Sensor("STD",     'Normal', 0, 0);

% asv = asv.Add_ParticleFilter(bathy);


% --- Diplay Initial Setup ------------------------------------------------
% fig = bathy.Plot_Map_UTM;
fig = ti.Map.PlotLayer_Contour("STD");

asv.Nav.PlotWaypoints(fig);

% fig = bathy.Plot_3DModel(0, 90);                                          % Plot Bathymetry Model
asv.PlotLocation(fig);


%% Run Simulation
speed = 3;              % ASV Speed                         [m/s]

sim_time = 0;           % Simulation time                   [s]
dt = 0.5;               % Simulation time step              [s]
tic
while toc < 1000 && ~asv.Goal_Achived(0.5) && ~asv.Nav.Done
    
    % ___ ASV Operations __________________________________________________
%     asv = asv.Move(goal, speed, dt);
    asv = asv.Navigate(speed, dt , 3);
    asv = asv.UpdatePF(dt, bathy.Depth(asv.State.XYZ.DR) );
    
    asv = asv.Add_SensorMeasurment("Entropy", ti.Map.Measure(asv.State.XYZ.GT(1:2), "Entropy"), sim_time);
    asv = asv.Add_SensorMeasurment("STD",     ti.Map.Measure(asv.State.XYZ.GT(1:2), "STD"),     sim_time);
    asv = asv.Add_SensorMeasurment("Depth",  -bathy.Depth( asv.State.XYZ.GT(1:2)),              sim_time);
    
    sim_time = sim_time +dt;                % time keeping
    
%     asv.PlotLocation(fig);
    
    if asv.Goal_Achived(3)
        asv.Nav.PF.PlotState(fig);
        asv.PlotPaths(fig);
    end
    
end

hold off

% asv.PlotPaths(fig);

% legend([p1,p2], "Start", "End")


%% Plotting
name_ = "LNH Low STD - Windy";
cov = asv.Nav.Cov;

t = zeros(size(cov));
for ii = 1: size(cov,2), t(ii) = trace(cov(ii).PF); end

h = asv.Sensors.Entropy.value;
s = asv.Sensors.STD.value;
d = asv.Sensors.Depth.value;

gt = asv.Paths.GT(1:asv.Paths.ind-1, 1:2);
pf = asv.Paths.PF(1:asv.Paths.ind-1, 1:2);

err = rmse(gt - pf);

figure('name', 'Covariance Trace')
subplot(4,1,1)
plot(t)
title(sprintf("Particle Filter Covariance Trace, %s Path",name_))

subplot(4,1,2)
plot(err)
title('Particle Filter Mean Squar Error')

subplot(4,1,3)
plot(s)
title('Terrain STD')

subplot(4,1,4)
plot(d)
title('Depth Profile')






