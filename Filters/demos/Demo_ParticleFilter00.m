close all
clear all
clc

ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(false);                                                      % Indicate whether new data should be selescted by the user

parmSweep = false;


%% Get Data
% Get files with file path GUI ----------------------------------------------------------------------------------------------
emfile = ui.GetFile('Ecomapper_Log','log');                                 % GUI to get the name and file path of a file
bathFile = ui.GetFile('Bathymetry_Map','mat');                          % GUI to get the name and file path of a file

if isempty(emfile) || isempty(bathFile)                                 % Check if the Data Input box was cnaceled
    return                                                                  % End script if there isn't an input for it to use
end

EM_Data = EcoMapperLog2Mat(emfile,'dvl filter');

            
            
if isempty(EM_Data)                                                         % Check that data is present
    disp("Empty Data Structure")
    disp("Doulbe check the directroy that you chose")
    return                                                                  % End script if there isn't any data to process
end

Bathymetry_Map = load(bathFile);                                            % This will load variables: Bathymetry, LON, LAT, Var
Bathymetry_Map = Bathymetry_Map.bathy_map;

[ emUTMs(:,1),  emUTMs(:,2), ~] = deg2utm( EM_Data.vehicle(:,1),    EM_Data.vehicle(:,2));     % Convert lat-lon to utms
[altUTMs(:,1), altUTMs(:,2), ~] = deg2utm( EM_Data.bathymetry(:,1), EM_Data.bathymetry(:,2));


%% Start Particle filter
pf = ParticleFilter00(10000, Bathymetry_Map);                                 % Instanciate the Particle filter with 1000 particles ane the Bathymetery Map

pf = pf.SetLocation( [emUTMs(1,1), emUTMs(1,2), EM_Data.attitude(1, 3)]);   % Set the initial vehilce location
pf = pf.SetNoise('sensor', 'Normal', 0, 0.9);                               % Set Sensor model as a probability distribuition ( "type of distribuition' , mu, sigma )
pf = pf.SetNoise('vehicle','Normal', 0, 0.2);                               % Set vehicle transiton model as a probability distribuition ( "type of distribuition' , mu, sigma )
pf = pf.SetNoise('compass','Normal', 0, 3);                                 % Set Comapss Noise

fig = Bathymetry_Map.Plot_3DModel(0, 90);                                   % Show Partilces over the Bathymetry Map
hold on

path = zeros(size(EM_Data.timeStamp,1),2);
path(1,:) = emUTMs(1,:);

for iter = 2: size(EM_Data.timeStamp,1)                                     % Run the partilce filter on the Data
       
    speed       = EM_Data.attitude(iter, 5);
    heading     = EM_Data.attitude(iter, 3);
    dT          = seconds( EM_Data.timeStamp(iter) - EM_Data.timeStamp(iter-1) );
    measurement = EM_Data.attitude(iter, 4);
    offSet      = altUTMs(iter,:) -  emUTMs(iter,:);
    
    pf = pf.Update(speed, heading, dT, measurement, offSet);
    
    path(iter,:) = pf.GetLocation;                                          % Get Vehicle location from Particle Filter
    
    plot(emUTMs(1:iter,1), emUTMs(1:iter,2), 'g')                           % Show Particles with Ecomapper Path at a straight on view
    plot(path(1:iter,1),path(1:iter,2), 'b')
    drawnow
    pause(0.25)
end




%% Plot Results

err = rms( emUTMs - path, 2 );

figure('name', 'Particle Filter Results', 'numbertitle', 'off')
subplot(1,2,1)
p1 = plot( emUTMs(:,1), emUTMs(:,2), 'g');
hold on
plot( emUTMs(end,1), emUTMs(end,2), '*g')
plot( emUTMs(end,1), emUTMs(end,2), 'og')

p2 = plot( path(:,1),   path(:,2), 'b');
plot( path(end,1), path(end,2), '*b')
plot( path(end,1), path(end,2), 'ob')
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
        
        pf = ParticleFilter(10000, Bathymetry_Map);                         % Instanciate the Particle filter with 1000 particles ane the Bathymetery Map
        
        pf = pf.SetLocation( [emUTMs(1,1), emUTMs(1,2)] );                  % Set the initial vehilce location
        pf = pf.SetNoise_Sensor( 'Normal', 0, params(jj) );                 % Set Sensor model as a probability distribuition ( "type of distribuition' , mu, sigma )
        pf = pf.SetNoise_Vehicle('Normal', 0, params(ii) );                 % Set vehicle transiton model as a probability distribuition ( "type of distribuition' , mu, sigma )
        
        path = zeros(s,2);
        path(1,:) = emUTMs(1,:);
        
        for iter = 2: size(s,1)                                             % Run the partilce filter on the Data
            
            speed       = EM_Data.attitude(itter, 5);                       % Deadreckoning values
            heading     = EM_Data.attitude(iter, 3) * 2 * pi / 360;
            dT          = seconds( EM_Data.timeStamp(iter) - EM_Data.timeStamp(iter-1) );
            measurement = EM_Data.attitude(iter, 4);
            
            offSet      = altUTMs(iter,:) -  emUTMs(iter,:);                % Off set between vehilce reading and bathymetery point due to vehicle attitude
            
            pf = pf.Update(speed, heading, dT, measurement, offSet);        % Up date particle filter with the new movememnt and measurment
            
            path(iter,:) = pf.GetLocation;                                  % Get Vehicle location from Particle Filter
            
        end
        
        err(jj) = rms(rms( emUTMs - path ));
        
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
pf = ParticleFilter(10000, Bathymetry_Map);                                  % Instanciate the Particle filter with 1000 particles ane the Bathymetery Map

pf = pf.SetLocation( [emUTMs(1,1), emUTMs(1,2)] );                          % Set the initial vehilce location
pf = pf.SetNoise_Sensor('Normal', 0, params(col));                                  % Set Sensor model as a probability distribuition ( "type of distribuition' , mu, sigma )
pf = pf.SetNoise_Vehicle('Normal', 0, params(row));                                 % Set vehicle transiton model as a probability distribuition ( "type of distribuition' , mu, sigma )

% pf.Plot                                                                     % Show Partilces over the Bathymetry Map

path = zeros(size(EM_Data.date,1),2);
path(1,:) = emUTMs(1,:);

for iter = 2: size(EM_Data.date,1)                                            % Run the partilce filter on the Data
       
    speed       = EM_Data.attitude(iter, 5);
    heading     = EM_Data.attitude(iter, 3) * 2 * pi / 360;
    dT          = seconds( EM_Data.date(iter) - EM_Data.date(iter-1) );
    measurement = EM_Data.attitude(iter, 4);
    offSet      = altUTMs(iter,:) -  emUTMs(iter,:);
    
    pf = pf.Update(speed, heading, dT, measurement, offSet);
    
    path(iter,:) = pf.GetLocation;                                            % Get Vehicle location from Particle Filter
    
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

plot( path(:,1),   path(:,2), 'b')
plot( path(end,1), path(end,2), '*b')
plot( path(end,1), path(end,2), 'ob')
hold off
drawnow
% export_fig CorrectedPath.png           

