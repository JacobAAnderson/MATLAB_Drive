close all
clear all
clc

format long g

ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(true);                                                      % Indicate whether new data should be selescted by the user

parmSweep = false;


%% Get Data
% Get files with file path GUI ----------------------------------------------------------------------------------------------
emfile = ui.GetFile('Ecomapper_Log','log');                                 % GUI to get the name and file path of a file
bathFile = ui.GetFile('Bathymetry_Map','mat');                              % GUI to get the name and file path of a file

if isempty(emfile) || isempty(bathFile)                                     % Check if the Data Input box was cnaceled
    return                                                                  % End script if there isn't an input for it to use
end

EM_Data = EcoMapperLog2Mat(emfile,'dvl filter');

            
if isempty(EM_Data)                                                         % Check that data is present
    disp("Empty Data Structure")
    disp("Doulbe check the directroy that you chose")
    return                                                                  % End script if there isn't any data to process
end

map = load(bathFile);                                                       % This will load variables: Bathymetry, LON, LAT, Var
map = map.bathy_map;


%% Set Up Particle Filter
[ emUTMs(:,1),  emUTMs(:,2), ~] = deg2utm( EM_Data.vehicle(:,1),    EM_Data.vehicle(:,2));     % Convert lat-lon to utms
[altUTMs(:,1), altUTMs(:,2), ~] = deg2utm( EM_Data.bathymetry(:,1), EM_Data.bathymetry(:,2));

pf = ParticleFilter('Demo', 10000);                                         % Instanciate the Particle filter with 1000 particles ane the Bathymetery Map
pf = pf.Add_Map(map);
pf = pf.SetLocation( [emUTMs(1,1), emUTMs(1,2), EM_Data.attitude(1, 3)]);   % Set the initial vehilce location

pf = pf.SetNoise('speed',  'Normal', 0, 0.75);                              % Set vehicle transiton model as a probability distribuition ( "type of distribuition' , mu, sigma )
pf = pf.SetNoise('compass','Normal', 0, 10);                                % Set Comapss Noise
pf = pf.SetNoise('alt',    'Normal', 0, 0.5);                               % Set Sensor model as a probability distribuition ( "type of distribuition' , mu, sigma )
pf = pf.SetNoise('acoms',  'Normal', 0, 0.5);                               % Set Sensor model as a probability distribuition ( "type of distribuition' , mu, sigma )


path_DR = zeros(size(EM_Data.timeStamp,1),2);
path_DR(1,:) = emUTMs(1,:);

path_PF = zeros(size(EM_Data.timeStamp,1),4);
path_PF(1,:) = [emUTMs(1,:), 0, 0];

x1   = [361972.04, 3701687.86];
cov1 = [ 1, 0.2; 0.2, 1];

x2   = [361855.48, 3701542.86];
cov2 = [ 1, -0.3; -0.3, 1];

x3   = [362016.84, 3701515.87];
cov3 = [ 1, 0; 0, 1];


%% Show Initial Configuration
fig = map.Plot_3DModel(0, 90);                                              % Show the Bathymetry Map
hold on

plot(x1(1), x1(2), '*b')
plot(x2(1), x2(2), '*b')
plot(x3(1), x3(2), '*b')


%% Run the route
for iter = 2: size(EM_Data.timeStamp,1)                                     % Run the partilce filter on the Data
       
    speed       = EM_Data.attitude(iter, 5);
    heading     = EM_Data.attitude(iter, 3);
    dT          = seconds( EM_Data.timeStamp(iter) - EM_Data.timeStamp(iter-1) );
    alt = EM_Data.attitude(iter, 4);

    pf = pf.Update(speed, dT, heading, alt);

    % ____ Range mesasurment updates ______________________________________
%     dist1 = vecnorm(emUTMs(iter,:) - x1) + random('Normal', 0, 0.5);
%     dist2 = vecnorm(emUTMs(iter,:) - x2) + random('Normal', 0, 0.5);
%     dist3 = vecnorm(emUTMs(iter,:) - x3) + random('Normal', 0, 0.5);
%     
%     acoms1 = {'acoms', dist1, x1, cov1 .* 0.5};
%     acoms2 = {'acoms', dist2, x2, cov2 .* 0.5};
%     acoms3 = {'acoms', dist3, x3, cov3 .* 0.5};
    
%     pf = pf.Update(speed, dT, heading, alt, acoms1, acoms2, acoms3);
    
    path_PF(iter,:) = pf.X';                                                % Get Vehicle location from Particle Filter
    path_DR(iter,:) = path_DR(iter-1,:) + dT * speed .* [sind(heading), cosd(heading)]; % Dead reckoning movement
    
    % ____ Plot Stuff _____________________________________________________
    plot(emUTMs(1:iter,1), emUTMs(1:iter,2),  'g')                          % Show Particles with Ecomapper Path at a straight on view
    plot(path_PF(1:iter,1),path_PF(1:iter,2), 'b')
    plot(path_DR(1:iter,1),path_DR(1:iter,2), 'r')
    drawnow
     
    if mod(iter, 100) == 0, pf.PlotState(fig, 0.5); hold on, end
    
end




%% Plot Results

% err = rms( emUTMs - path(:,1:2), 2 );

err = rmse( emUTMs - path_PF(:,1:2));

figure('name', 'Particle Filter Results', 'numbertitle', 'off')
subplot(1,2,1)
p1 = plot( emUTMs(:,1), emUTMs(:,2), 'g');
hold on
plot( emUTMs(end,1), emUTMs(end,2), '*g')
plot( emUTMs(end,1), emUTMs(end,2), 'og')

p2 = plot( path_PF(:,1),   path_PF(:,2), 'b');
plot( path_PF(end,1), path_PF(end,2), '*b')
plot( path_PF(end,1), path_PF(end,2), 'ob')
hold off

xlabel('Easting')
ylabel("Northing")
title('GPS Path vs Partilce Filter Path')
legend([p1,p2],{'GPS', 'PF'})


subplot(1,2,2)
plot(EM_Data.timeStamp, err, 'r')
xlabel('Date-Time')
ylabel('RMS Error [m]')
title('Particle Filter RMS Error')
drawnow
% export_fig CorrectedPath.png                                              % Save plot as figure 


%% Built in Parameter Sweep ???
% params = 0:0.1:2; 
% newParams = pf.parameterSweep( params, EM_Data, true);


%% Parameter sweep
if ~parmSweep
    return
end


disp(" Doing Parameter sweep ")
tic
params = 0:0.1:2;                                                           % Range of Vehicle uncertainty

s = size(EM_Data.timeStamp,1);

parfor ii = 1: numel(params)                                                % Iterate throught 
    
    err = zeros(size(params));                                              % Prealocate space for error vector
    
    for jj =1: numel(params)
        
        pf = ParticleFilter(10000, map);                         % Instanciate the Particle filter with 1000 particles ane the Bathymetery Map
        
        pf = pf.SetLocation( [emUTMs(1,1), emUTMs(1,2)] );                  % Set the initial vehilce location
        pf = pf.SetNoise_Sensor( 'Normal', 0, params(jj) );                 % Set Sensor model as a probability distribuition ( "type of distribuition' , mu, sigma )
        pf = pf.SetNoise_Vehicle('Normal', 0, params(ii) );                 % Set vehicle transiton model as a probability distribuition ( "type of distribuition' , mu, sigma )
        
        path_PF = zeros(s,2);
        path_PF(1,:) = emUTMs(1,:);
        
        for iter = 2: size(s,1)                                             % Run the partilce filter on the Data
            
            speed       = EM_Data.attitude(itter, 5);                       % Deadreckoning values
            heading     = EM_Data.attitude(iter, 3) * 2 * pi / 360;
            dT          = seconds( EM_Data.timeStamp(iter) - EM_Data.timeStamp(iter-1) );
            measurement = EM_Data.attitude(iter, 4);
            
            offSet      = altUTMs(iter,:) -  emUTMs(iter,:);                % Off set between vehilce reading and bathymetery point due to vehicle attitude
            
            pf = pf.Update(speed, heading, dT, measurement, offSet);        % Up date particle filter with the new movememnt and measurment
            
            path_PF(iter,:) = pf.GetLocation;                                  % Get Vehicle location from Particle Filter
            
        end
        
        err(jj) = rms(rms( emUTMs - path_PF ));
        
    end
    
    Err(ii,:) = err;
    
end
toc
disp("Done!!")

%% Show Resutls from Parameter sweep
close all
[row, col] = find(Err == min(Err(:)));
[X, Y] = meshgrid(params, params);

figure('Name', "Parameter Sweep Results")
surf(X,Y,Err);
colorbar
hold on
p = plot3(params(col), params(row), Err(row,col), '*r');
hold off

xlabel("Sensor Noise")
ylabel("Vehicle Noise")
zlabel("RMS Error")
legend(p,"Mininum RMS Err",'location','north')
title(sprintf("RMS Error, Min at --> X: %f, Y: %f", params(col), params(row)))
view(0,90)



%% Run Particle filter again with resutls from parameter sweep
pf = ParticleFilter(10000, map);                                  % Instanciate the Particle filter with 1000 particles ane the Bathymetery Map

pf = pf.SetLocation( [emUTMs(1,1), emUTMs(1,2)] );                          % Set the initial vehilce location
pf = pf.SetNoise_Sensor('Normal', 0, params(col));                                  % Set Sensor model as a probability distribuition ( "type of distribuition' , mu, sigma )
pf = pf.SetNoise_Vehicle('Normal', 0, params(row));                                 % Set vehicle transiton model as a probability distribuition ( "type of distribuition' , mu, sigma )

% pf.Plot                                                                     % Show Partilces over the Bathymetry Map

path_PF = zeros(size(EM_Data.date,1),2);
path_PF(1,:) = emUTMs(1,:);

for iter = 2: size(EM_Data.date,1)                                            % Run the partilce filter on the Data
       
    speed       = EM_Data.attitude(iter, 5);
    heading     = EM_Data.attitude(iter, 3) * 2 * pi / 360;
    dT          = seconds( EM_Data.date(iter) - EM_Data.date(iter-1) );
    measurement = EM_Data.attitude(iter, 4);
    offSet      = altUTMs(iter,:) -  emUTMs(iter,:);
    
    pf = pf.Update(speed, heading, dT, measurement, offSet);
    
    path_PF(iter,:) = pf.GetLocation;                                            % Get Vehicle location from Particle Filter
    
%     close all
%     pf.Plot( [emUTMs(1:ii,1), emUTMs(1:ii,2)], [0,90] )                    % Show Particles with Ecomapper Path at a straight on view
%     pause(0.25)
end
 
% Plot Results
figure
plot( emUTMs(:,1), emUTMs(:,2), 'r')
hold on
plot( emUTMs(end,1), emUTMs(end,2), '*r')
plot( emUTMs(end,1), emUTMs(end,2), 'or')

plot( path_PF(:,1),   path_PF(:,2), 'b')
plot( path_PF(end,1), path_PF(end,2), '*b')
plot( path_PF(end,1), path_PF(end,2), 'ob')
hold off
drawnow
% export_fig CorrectedPath.png           

