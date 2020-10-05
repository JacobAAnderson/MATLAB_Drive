%% Riptide Dec-TBN
% Jacob Anderson
% Robotic Decition Making Laboratory (RDML)
% Oregon State University
% Corallis OR
% August 25, 2020

close all
clear all
clc


%% Get Data
ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(false);                                                      % Indicate whether new data should be selescted by the user

[RT, geotiff, bathy] = GetRiptides(ui, 2, {'210','216'});                   % Choose data files with gui and return data structure

dtfs = DateTime_Funs;                                                       % Date Time Functions


%% Prep Data

% fig1 = bathy.Plot_3DModel(0, 90);                                           % Show the Bathymetry Map

% ---- Initialization Sequence -----

% missions = [3,2];         % Test Lake Mission Set
missions = [1,3];           % Fosters Lake Mission Set

for v = 1: numel(RT)
    RT(v) = RT(v).Select_Mission( missions(v) );                            % Choose Mission
    RT(v) = RT(v).Model_VehicleSpeed;                                       % Calculate corrected speed
    RT(v) = RT(v).Get_VehiclePaths;                                         % Evaluate GT and DR paths
    RT(v) = RT(v).Add_ParticleFilter(bathy);                                % Prep Particle filters
    
%     fig1 = RT(v).pf.PlotState(fig1);                                        % Display Particle filter state
    
end

clear missions v startime

RT_ = RT;                                                                   % Keep a clean coppy of the vehicle data


% ---- Model the latency in acoustic communications ---
[mu, sig] = Model_AcommsLatency( RT );


% --- Show Acoustic Communications ---
[~] = Plot_AcousticCommunications(RT, geotiff);

% --- Show vehicle paths ---
% fig1 = geotiff.Show;
% 
% for rt = RT
%     fig1 = rt.Plot_Paths(fig1);
%     fig1 = rt.Plot_Course(fig1);
%     fig1 = rt.Plot_Acomms(fig1);
% end

% fig1 = RT(2).Plot_Path(fig1, 'DR', 'y');
% fig1 = RT(2).Plot_Path(fig1, 'GT', 'w');
% fig1 = RT(2).Plot_Course(fig1);









return

%% Do TBN

% -- Fresh Start --
close all   
clear RT time v mem in speed heading dt measurement s rx_msg rx_time tx_time owtof dist x sig acoms rt idx fig1 fig2
RT = RT_;

% ---- Get Timeline -----
timeLine = [RT(1).data.gpsDate;
            RT(1).acomms.gpsDate;
            RT(2).data.gpsDate;
            RT(2).acomms.gpsDate];

timeLine = sort(unique(timeLine));

count = [0,0];

% ---- Run through data logs doing TBN -----
disp('Performing TBN')
for time = timeLine'                                                        % Run the timeline
    for v = 1:numel(RT)                                                     % Switch between vehicles
        for mem = {'data','acomms'}                                         % Switch between data sources
            for in = find(dtfs.Ismember(RT(v).(mem{1}).gpsDate, time))'          % Find which data it is time for
                
                % ---- Process Navigation data ----
                if strcmp(mem{1},'data')
                    
                    % --- Particle Filter Update ---
                    speed       = RT(v).filteredData.speed(in);              % Use the updated speed 
%                     speed       = RT(v).data.attitude(in,5);              % Use raw speed
                    heading     = RT(v).data.attitude(in,3);
                    dt          = RT(v).data.timeStep(in);
                    measurement = RT(v).data.bathymetry(in,3);
                    
                    RT(v).pf = RT(v).pf.Update(speed, dt, heading, measurement);
                   
                    count(v) = count(v) + 1;
                    
                end
                
                
                
                % ---- Process Recived Acomms messages ----
                if strcmp(mem{1},'acomms') && ~isnat(RT(v).acomms.recived_msg(in))
                    
                    if v == 1, s=2;                                         % Get index of the vehicle that sent the message
                    else, s = 1;
                    end
                    
                    rx_msg  = RT(v).acomms.recived_msg(in);                 % Get the recived message
                    
                    idx = dtfs.Ismember(RT(s).acomms.sent_msg, rx_msg);          % Match it to the message sent from the other vehicle
                    
                    if ~any(idx), continue, end                             % Continue if there is not a matching sent message
                    
                    tx_time = RT(s).acomms.gpsDate(idx);
                    rx_time = RT(v).acomms.gpsDate(in);
                    owtof   = seconds(rx_time - tx_time);                   % Compute one way time of flight
                    
                    dist  = owtof * 1475;
                    x     = RT(s).pf.X(1:2)';
                    sig   = RT(s).pf.Cov(1:2,1:2);
                    acoms = {'Acoms', dist, x, sig};
                    
                    RT(v).pf = RT(v).pf.Acoms(acoms);
                    
                    
                end
                
                
            end
        end
    end
end

disp(count)

% ____ Post Processing ____
for v = 1: numel(RT)
    RT(v) = RT(v).Get_PFpath;                                               % Copy over PF path to RT and conver to lat_lon 
end


% % ____ Display results ____
% fig1 = bathy.Plot_3DModel(0, 90);                                           % Show the Bathymetry Map
% fig2 = geotiff.Show;                                                        % Show Geotiff
% 
% for rt = RT
%     
%     fig1 = rt.pf.PlotPath(fig1);                                            % Plot utm path on bathymetry map
%     fig2 = rt.Plot_Paths(fig2);                                             % Plot lat-lon path on geotiff
%     
% end


clear time v mem in speed heading dt measurement s rx_msg rx_time tx_time owtof dist x sig acoms rt idx



%% ---- Calculate path errors and generate figures -----

paths = {'DR_corrected', 'PF'};
types = {'Error_210', 'Error_216', 'JointErr', 'Time'};

rAd = ExtractPathData(RT, paths, types);

for type = types(1:3), rAd.PlotData( paths, type{1}, 'Time'), end

fprintf('\n\n\n')
fprintf("Total Joint Error for Dear Reckoning:  %.2f [m s]\n", rAd.Stats.DR_corrected.JointErr.Area_Ave)
fprintf("Total Joint Error for Particle Filter: %.2f [m s]\n", rAd.Stats.PF.JointErr.Area_Ave)
fprintf('\n\n\n')

clear paths types type


%% Display Specific Graphs
fig1 = RT(1).Plot_Paths(geotiff);                                           % Plot lat-lon path on geotiff
fig2 = RT(2).Plot_Paths(geotiff);                                           % Plot lat-lon path on geotiff



%% Functions

function rAd = ExtractPathData(RT, paths, types)

num   = min(size(RT(1).path.GT.lat_lon,1), size(RT(2).path.GT.lat_lon,1));

rAd = RandomAssData([num,1], paths, types);                                 % Instanciate RandomAssData

for sim = paths
    
    err1 = abs( RT(1).path.(sim{1}).utm(1:num,:) - RT(1).path.GT.utm(1:num,:) );
    err2 = abs( RT(2).path.(sim{1}).utm(1:num,:) - RT(2).path.GT.utm(1:num,:) );
    jointErr = err1 + err2;
    
    for time = 1:num
        rAd = rAd.Add_Result(1, sim{1}, err1(time), err2(time), jointErr(time), time); % Add data to the respective set
    end
end

rAd = rAd.Eval_Stats('Time');

end

