%% Riptide Dec-TBN
% Jacob Anderson
% Robotic Decition Making Laboratory (RDML)
% Oregon State University
% Corallis OR
% August 25, 2020

close all
clear all
clc

dtfs   = DateTime_Funs;                                                     % Date Time Functions
rtAcms = Riptide_Acomms;                                                    % Joint Riptide functions
sos    = 1475;                                                              % Speed of sound in water [m/s];



%% Get Data
ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(false);                                                     % Indicate whether new data should be selescted by the user

filters = {'GPS_fix', 0};                                                   % Data filter {"Data_Field", value to be filterred out}
%            'ALT_ALTITUDE', 0};        

[RT, geotiff, bathy] = GetRiptides(ui, 2, {'Dory','Nemo'}, filters);        % Choose data files with gui and return data structure



% Sort Riptide_Acomms
[~, ~, owtof] = rtAcms.Model_AcommsLatency( RT );                            % Model the latency in acoustic communications and get one way time of flight

for v = 1:numel(RT)
    RT(v) = RT(v).Add_OWTOF(owtof);
    RT(v) = RT(v).Get_Manifest("skipidel", "alt & acomms");
    
    RT(v).commsPlanner = CommunicationPlanning(RT(v).name);
end



%% ____ Select Mission ____
close all
clc

% missions = [ 1, 1];         % Fosters Lake Mission Set: CurcumNav1, CircumNav2
% missions = [ 2, 2];         % Fosters Lake Mission Set: CurcumNav1, CircumNav2        --> Good Acoms
% missions = [ 3, 3];         % Fosters Lake Mission Set: Dimond1_A, Dimond1_B          --> Decent Acomms
% missions = [ 5, 7];         % Fosters Lake Mission Set: test1, test2
% missions = [ 7, 8];         % Fosters Lake Mission Set: ZZ1, ZZ2
% missions = [ 8, 9];         % Fosters Lake Mission Set: ZZ1, ZZ2
% missions = [ 9,10];         % Fosters Lake Mission Set: ZZ3_A, ZZ3_B                  --> Decent Acoms
% missions = [10,11];         % Fosters Lake Mission Set: ZZ3_A, ZZ3_B                  --> Decent Acomms
% missions = [11,12];         % Fosters Lake Mission Set:CircumNav2_B, CircumNave2_A    --> Decent Acomms
 missions = [14,15];         % Fosters Lake Mission Set:ZZ1, ZZ2                         --> Decent Acomms
% missions = [15,16];         % Fosters Lake Mission Set:Dimond1_A, Dimond1_B  (DEC-TBN Does better)

fig1 = geotiff.Show;

for v = 1: numel(RT)
    
    RT(v) = RT(v).Select_Mission( missions(v) );                            % Choose Mission
    RT(v) = RT(v).Get_VehiclePaths;                                         % Evaluate GT and DR paths for speed and compass modeling
    RT(v) = RT(v).Model_VehicleSpeed;                                       % Calculate corrected speed from GPS
    RT(v) = RT(v).Model_Compass('poly2');                                   % Model / Calibrate Compass data
    RT(v) = RT(v).Model_Altimiter(20);                                      % Filter Altimiter data and build normal distribution
    RT(v) = RT(v).Get_VehiclePaths;                                         % Evaluate GT and DR paths
    RT(v) = RT(v).Add_TBN(bathy);                                           % Prep Particle filters
        
    RT(v).Disp_Manifest;                                                    % Display Vehicle Manifest
    
    %___ Plot Stuff ____
    fig1 = RT(v).Plot_Paths(fig1);                                          % Show paths on a geotiff
%     [~] = RT(v).Plot_Altitude;                                              % Show Altimeter profile
%     [~] = RT(v).Plot_Alt_Geotiff(geotiff);                                  % Show Altimeter location on map

end


% --- Get Acoustic modem model for the selected mission ---
[mu, sig] = rtAcms.Mission_Latency(RT);                                      

% - Set Acomms Noise Distributions -
for v = 1:numel(RT)
    m = seconds(mu(v));
    s = seconds(sig(v))*sos;
    RT(v).tbn = RT(v).tbn.SetNoise('acoms',  'Normal', m, s);                 % Set acoustic modem error Distributions
end


[~] = rtAcms.Plot_AcousticCommunications(RT, geotiff);                      % Plot Acoustic Communications

clear missions v owtof filters sig m s

RT_ = RT;                                                                   % Keep a clean coppy of the vehicle data



%% ---- Get Timeline -----
disp("--- Getting Time Line ---")
tic
timeLine = [RT(1).data.gpsDate;
            RT(1).acomms.gpsDate;
            RT(2).data.gpsDate;
            RT(2).acomms.gpsDate];

timeLine = sort(dtfs.Unique(timeLine));
fprintf("Time Line finished\nElapsed time: %.1f [sec]\n\n", toc)


%% ____ Do TBN ____

% -- Fresh Start --
% clc
close all
clear RT time v mem in speed heading dt measurement s rx_msg rx_time tx_time owtof dist x sig acoms rt idx fig1 fig2
RT = RT_;


count = [0,0];

% ---- Run through data logs doing TBN -----
disp('==== Performing TBN ====')
for time = timeLine'                                                        % Run the timeline
    for v = 1:numel(RT)                                                     % Switch between vehicles
        for mem = {'data'}%,'acomms'}                                       % Switch between data sources
            for in = find(dtfs.Ismember(RT(v).(mem{1}).gpsDate, time))'     % Find which data it is time for
                
                % ---- Process Navigation data ----
                if strcmp(mem{1},'data')
                    
                    % --- Particle Filter Update ---
                    speed   =  RT(v).filteredData.speed(in);                 % Use the updated speed 
                    heading =  RT(v).data.attitude(in,3);
                    dt      =  RT(v).data.timeStep(in);
                    alt     = -RT(v).filteredData.altitude(in);
                    
%                     fprintf("[%s]  TBN - Update, speed: %f, heading: %f, dt: %f", RT(v).name, speed, heading, dt) 
                    
%                     disp(" ")
%                     RT(v).tbn = RT(v).tbn.Update(speed, dt, heading);
                    
                    if alt == 0
%                         fprintf('\n')
                        RT(v).tbn = RT(v).tbn.Update(speed, dt, heading);
                    else
%                         fprintf(', alt: %f\n', alt)
                        RT(v).tbn = RT(v).tbn.Update(speed, dt, heading, alt);
                    end
%                     count(v) = count(v) + 1;
                    
                end
                
                
                
                % ---- Process Recived Acomms messages ----
                if false strcmp(mem{1},'acomms') && ~isnat(RT(v).acomms.recived_msg(in))
                    
                    if v == 1, s=2;                                         % Get index of the vehicle that sent the message
                    else, s = 1;
                    end
                    
                    owtof = RT(v).acomms.tof(in) - mu(v);                   % One way time of flight
                    dist  = seconds(owtof) * sos;                           % Compute Distance
                    x     = RT(s).tbn.X(1:2)';
                    sig   = RT(s).tbn.cov(1:2,1:2);
%                     acoms = {'Acoms', dist, x, sig};
                    
                    fprintf("[%s]  TBN - Acomms, owtof: %1.8s, dist: %5.2f [m]\n",RT(v).name, owtof, dist) 
                    
                    RT(v).tbn = RT(v).tbn.Acoms(dist, x, sig);
                    count(v) = count(v) + 1;
                    
                end
                
                
            end
        end
    end
end

disp(count)

% ____ Post Processing ____
for v = 1: numel(RT)
    RT(v) = RT(v).Get_TBNpath;                                               % Copy over TBN path to RT and conver to lat_lon 
end


% % ____ Display results ____
% fig1 = bathy.Plot_3DModel(0, 90);                                           % Show the Bathymetry Map
% fig2 = geotiff.Show;                                                        % Show Geotiff
% 
% for rt = RT
%     
%     fig1 = rt.tbn.PlotPath(fig1);                                            % Plot utm path on bathymetry map
%     fig2 = rt.Plot_Paths(fig2);                                             % Plot lat-lon path on geotiff
%     
% end


clear time v mem in speed heading dt measurement s rx_msg rx_time tx_time owtof dist x sig acoms rt idx



% ---- Calculate path errors and generate figures -----

paths = {'DR_corrected', 'tbn'};
types = {'Error_210', 'Error_216', 'JointErr', 'Time'};

rAd = ExtractPathData(RT, paths, types);


% _____ Display Results __________________________________________________
for type = types(1:3), rAd.PlotData( paths, type{1}, 'Time'), end

fprintf('\n\n\n')
fprintf("Total Joint Error for Dead Reckoning:  %.2f [m s]\n", rAd.Stats.DR_corrected.JointErr.Area_Ave)
fprintf("Total Joint Error for Particle Filter: %.2f [m s]\n", rAd.Stats.tbn.JointErr.Area_Ave)
fprintf('\n\n\n')


% ____ Display Specific Graphs ___________________________________________
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

