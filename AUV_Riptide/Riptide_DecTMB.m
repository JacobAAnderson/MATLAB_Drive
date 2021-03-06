%% Riptide Dec-TBN
% Jacob Anderson
% Robotic Decition Making Laboratory (RDML)
% Oregon State University
% Corallis OR
% August 25, 2020

close all
clear all
clc

ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(false);                                                     % Indicate whether new data should be selescted by the user


dtfs   = DateTime_Funs;                                                     % Dathlkhklhe Time Functions
rtAcms = Riptide_Acomms;                                                    % Joint Riptide functions
sos    = 1475;                                                              % Speed of sound in water [m/s];

timeLine_mission = [];                                                      % Keep track of which time line we have


%% Get Data
filters = {'GPS_fix', 0};                                                   % Data filter {"Data_Field", value to be filterred out}
%            'ALT_ALTITUDE', 0};        

[RT, geotiff, bathy] = GetRiptides(ui, 2, {'Dory','Nemo'}, filters);        % Choose data files with gui and return data structure


% Sort Riptide_Acomms
close all
clc
[~, ~, owtof] = rtAcms.Model_AcommsLatency( RT , ["CircumNav2_B", "CircumNav2_A", "idel"] );                           % Model the latency in acoustic communications and get one way time of flight

for v = 1:numel(RT)
    RT(v) = RT(v).Add_OWTOF(owtof);                                         % Add one way time of fligth 
    RT(v) = RT(v).Get_Manifest("skipidel", "alt & acomms");
    
    RT(v).commsPlanner = CommunicationPlanning(RT(v).name);
end

disp('Saving Riptide Data')
save("/Users/jake/_OutPuts/RT_FostersLake.mat", 'RT', 'geotiff', 'bathy')

clear owtof filters



%% Select Mission
close all
clc

load('/Users/jake/_OutPuts/RT_FostersLake.mat')                             % Start with a clean copy of the data

% ------ Missions   2020-10-27   ------------------------------------------
% missions = [1,1];  name_ = "Mission 0";           %    Fosters Lake Mission Set:CircumNav2_B, CircumNave2_A   --> Decent Acomms
% missions = [3,3];  name_ = "Mission 1";           % *! Fosters Lake Mission Set:Dimond1_A, Dimond1_B          --> Good TBN
missions = [4,4];  name_ = "Mission 2";           %    Fosters Lake Mission Set:ZZ1, ZZ2                      --> Decent Acomms
% missions = [5,5];  name_ = "Mission 3";           % *! Fosters Lake Mission Set:Dimond1_A, Dimond1_B          --> (DEC-TBN Does better)
% missions = [6,6];  name_ = "Mission 4";           % *!! Fosters Lake Mission Set:ZZ4_B, ZZ4_A  
 

% ----- Setup -------------------------------------------------------------
% fig1 = geotiff.Show;                                                      % Show geotiff for plotting paths

for v = 1: numel(RT)
    
    RT(v) = RT(v).VehicleModeling('poly2', bathy, missions(v) );            % Do vehicle modeling that excludes data from the selected mission
    
    % -- Select mission --
    RT(v) = RT(v).Select_Mission( missions(v) );                            % Choose Mission
    RT(v) = RT(v).Get_VehiclePaths;                                         % Evaluate GT and DR paths for speed and compass modeling
    
    % -- Vehicle Modeling --
    RT(v) = RT(v).Model_VehicleSpeed;                                       % Calculate corrected speed from GPS
    RT(v) = RT(v).Model_Compass('poly2');                                   % Model / Calibrate Compass data
    RT(v) = RT(v).Get_VehiclePaths;                                         % Evaluate GT and Corrected DR paths from the speed and compass models
    
    % -- Terrain Based Navigation stuff --
    RT(v) = RT(v).Add_TBN(bathy, 1);                                        % Prep Particle filters with std inflation ratio
    RT(v) = RT(v).Model_Altimiter(false);                                   % Filter Altimiter data and build normal distribution
   
    % -- Generic DEC-TBN Policies ---
    RT(v) = RT(v).Add_CommsPloicy( 'TBN', false);                           % Add TBN policy with all false
    RT(v) = RT(v).Add_CommsPloicy( 'Full', true);                           % Add Full DEC-TBP policy with all true
    
    % -- Plot/display Stuff --
%    RT(v).Disp_Manifest;                                                    % Display Vehicle Manifest
%    RT(v).Plot_Course(geotiff);
%     fig1 = RT(v).Plot_Paths(fig1);                                        % Show paths on a geotiff
%     RT(v).Plot_AltimeterProfile;                                          % Show Altimeter profile
%     RT(v).Plot_Altimeter(geotiff);                                        % Show Altimeter location on map

end

% RT(v) = RT(v).CorrectCourse(geotiff);                                       % Correct Nemo's course for the comms planning


[mu, sig] = rtAcms.Mission_Latency(RT);                                     % Get Acoustic modem model for the selected mission


for v = 1:numel(RT)                                                         % Set Acomms Noise Distributions 
    m = seconds(mu(v));
    s = seconds(sig(v))*sos;
    RT(v).tbn = RT(v).tbn.SetNoise('acoms',  'Normal', m, s * 1.5);         % Set acoustic modem error Distributions
end

% rtAcms.Plot_AcousticCommunications(RT, geotiff);                          % Plot Acoustic Communications

clear v m s sig




% ----- Get Timeline ------------------------------------------------------
if ~isequal( timeLine_mission, missions)
    disp("--- Getting Time Line ---")
    tic

    timeLine_mission = missions;
    timeLine = [RT(1).data.gpsDate;
                RT(1).acomms.gpsDate;
                RT(2).data.gpsDate;
                RT(2).acomms.gpsDate];
    
    timeLine = sort(dtfs.Unique(timeLine));
    fprintf("Time Line finished\nElapsed time: %.1f [sec]\n\n", toc)
end




% _____ Communication planning ____________________________________________
% clc
% close all


% ----- Tune Particle Filter for planning ---------------------------------
for prop = { 'altimeter', 'compass' 'speed';           % --> Type of noise
             'Normal',    'Normal', 'Normal';          % --> Default distribion type
              0,           0,        0;                % --> Default mu
              2.0,         25,       1.5;              % --> Sigma for Dory
              2.0,         25,       1.5}              % --> Sigma for Nemo
    

        RT(1).tbn = RT(1).tbn.SetNoise(prop{1}, prop{2}, prop{3}, prop{4});
        RT(2).tbn = RT(2).tbn.SetNoise(prop{1}, prop{2}, prop{3}, prop{5});

end



for ii = [1: numel(RT); numel(RT):-1:1], v = ii(1); s = ii(2);

    RT(v) = RT(v).Add_CommsPlanner;                                                     % Instanciate communications planner
    
    RT(v).commsPlanner = RT(v).commsPlanner.Add_Vehicle(RT(v).VehicleInfo('self'));     % Load vehicle info into the planner
    RT(v).commsPlanner = RT(v).commsPlanner.Add_Vehicle(RT(s).VehicleInfo);
   
end


cov_max = 1.25; 
threshold = 25;


dt = 30;                % Time step in seconds
dt_ = seconds(dt);      % time Step as a duration type

% _____ Plan Comms Policy _________________________________________________
% for ii = [1: numel(RT); numel(RT): -1: 1], v = ii(1); s = ii(2);
%     
%     fprintf('\nComms Planning for %s', RT(v).name)
%     
%     start = RT(v).data.timeStamp(1);
%     stop  = RT(v).data.timeStamp(end);
%     
%     sz = floor( (stop - start)/ dt_ );
%     
%     policy = [false; RT(v).commsPlanner.Plan2(sz, dt, cov_max, threshold, false)];        % Generate the policy and append a false to the front since it the planing starts at t = 2
%     
%     policyTimes = start: dt_: stop;
%     
%     RT(s) = RT(s).Add_CommsPloicy('Plan', policyTimes(policy), dt);         % Tell the other vehicle which communications it can use
%     
%     fprintf("Size: %d\n", numel(policy))
%     
% end

% _____ Generate Random Comms Policy ______________________________________

rando{6} = [];

for p = {0.8,    0.6,   0.4,    0.2,   0.1,   0.05;
        'eight' 'six', 'four', 'two', 'one', 'half';
         1,      2,     3,      4,     5,     6 }
    
    
    for ii = [1: numel(RT); numel(RT): -1: 1], v = ii(1); s = ii(2);
        
        fprintf('\n%f Random Comms for %s\n', p{1}, RT(v).name)
        
        start = RT(v).data.timeStamp(1);
        stop  = RT(v).data.timeStamp(end);
        
        sz = floor( (stop - start)/ dt_ );
        
        policy = false(sz+1,1);
        
        comms = randi(sz+1, floor(sz*p{1}), 1 );        % Generate the policy and append a false to the front since it the planing starts at t = 2
        
        policy(comms) = true;
        
        policyTimes = start: dt_: stop;
        
        RT(s) = RT(s).Add_CommsPloicy( sprintf('Random_%s', p{2}), policyTimes(policy), dt);         % Tell the other vehicle which communications it can use
        
        fprintf("Size: %d\t# of comms: %d\n", numel(policy), sum(policy) )
        
    end
    
    rando{p{3}} = sprintf('Random_%s', p{2});
    
end

clear ii v s dt dt_ policy policyTimes sz stop start


% _____ Do TBN ____________________________________________________________
%clc
%close all
 
num_trials = 1;

% --- Instanciate rAd to track data ---
% paths = [{'DR_Corrected', 'TBN', 'Full', 'Plan'}, rando];                  % Set up rAd with DR corrected data
paths = rando;
types = {'Dory_Error', 'Nemo_Error', 'Joint_Error', 'Time'};

[rAd, gtPaths] = GetRAD(RT, paths, types, num_trials, timeLine);            % Create rAd with DR Correctd data entered



% ----- Tune Particle Filter for Localizing -------------------------------
for prop = { 'altimeter', 'compass' 'speed';           % --> Type of noise
             'Normal',    'Normal', 'Normal';          % --> Default distribion type
              0,           0,        0;                % --> Default mu
              2.5,         30,       2.0;              % --> Sigma for Dory
              2.5,         30,       2.0}              % --> Sigma for Nemo
    

        RT(1).tbn = RT(1).tbn.SetNoise(prop{1}, prop{2}, prop{3}, prop{4});
        RT(2).tbn = RT(2).tbn.SetNoise(prop{1}, prop{2}, prop{3}, prop{5});

end

clear prop


% --- TBN Simulation ---
% for dec_tbn = {'TBN', 'Full', 'Plan', 'Random'}                             % Switch between TBN and DEC-TBN
for dec_tbn = paths                                                         % Switch between TBN and DEC-TBN
    
    for iter = 1: num_trials
        
        [tbnPaths, sim_times, ~] = Riptide_TBN_Sim(RT, timeLine, sos, mu, dec_tbn{1});    % Do TBN / DEC-TBN simulation
        
        rAd = TallyPathData(rAd, iter, dec_tbn{1}, gtPaths, tbnPaths, sim_times, timeLine);     % Calculate path errors and add to rAd
          
    end
    
    
    for v = 1: numel(RT)
        
        RT(v) = RT(v).Add_Path(tbnPaths(:,:,v), dec_tbn{1}, 'utm');
        
    end
    
    
end

clear dec_tbn iter tbnPaths sim_times


rAd = rAd.Eval_Stats('Time');                                               % Determine average error and total error


% ____ Display Graphs _____________________________________________________
close all
open('/Users/jake/OSU/RDML/Multi_Robot_Localization/Riptide/FostersLake/Figures/2020-11-21_Mission2_jointError.fig')

hold on

DisplyRAD(rAd, types{3}, paths)                                                % Display Result Stats

fprintf("Duration %.2f [sec]\n\n", seconds(timeLine(end) - timeLine(1)))


% for v = 1:numel(RT)
%     
%      RT(v).Plot_Paths(geotiff, [{'GT'}, paths], 'grbcm', name_);            % Plot lat-lon path on geotiff
%    
% end





%% ========= Functions ====================================================


%% Create rAd and out put GT paths for later err calculations
function [rAd, gtPaths] = GetRAD(RT, paths, types, num_trials, timeLine)

num = max( arrayfun( @(x) size( x.path.GT.utm, 1), RT));
rAd = RandomAssData([num,num_trials], paths, types);

gtPaths = zeros(num,2,numel(RT));

for v = 1: numel(RT)
    
    sz = size( RT(v).path.GT.utm, 1);
    gtPaths(1:sz,:,v) = RT(v).path.GT.utm;
    
end


drPaths = zeros(num,2,numel(RT));
drTimes = NaT(num, numel(RT));

for v = 1: numel(RT)
    
    sz = size( RT(v).path.GT.utm, 1);
    drPaths(1:sz,:,v) = RT(v).path.DR_Corrected.utm;
    
    drTimes(1:sz,v) = RT(v).data.gpsDate;
    
    
end

rAd = TallyPathData(rAd, 1, paths{1}, gtPaths, drPaths, drTimes, timeLine);


end




%% Add Path data to rAd 
function rAd = TallyPathData(rAd, iter, path, gtPaths, tbnPaths, times, t)


s = size(gtPaths,3);

err = NaN(size(t,1), s + 1 );



for ii = 1: size(t,1)
    
    for v = 1:s
        
        idx = find( ismember(times(:,v), t(ii)) );
        
        if any(idx)
            
            err(ii,v) = vecnorm( gtPaths(idx(1),1:2, v) - tbnPaths(idx(1),1:2, v), 2, 2);
            
            
        end
        
    end
    
end


for v = 1: s
    for ii = find(isnan(err(:,v)))'
        
        if ii == 1, continue, end
        
        err(ii,v) = err(ii-1,v);
        
    end
end


err(:,end) = nansum(err, 2);

for ii = 1:size(err,1)
    rAd = rAd.Add_Result(iter, path, err(ii,1), err(ii,2), err(ii,3), t(ii) ); % Add data to the respective set
end


rAd = rAd.ResetIndex;


end



%% Display rAd stats
function DisplyRAD(rAd, types, paths)

% for t = types(1:3), rAd.PlotData( paths, t{1}, 'Time'), end
% for t = types(1:3), rAd.PlotStat( paths, 'Ave', t{1}, 'Time'), end
rAd.PlotStat( paths, 'Ave', types, 'Time')

fprintf('\n\n\n')

for p = paths

    fprintf("Total Joint Error for %s:  %.2f [m s]\n", p{1}, rAd.Stats.(p{1}).Joint_Error.Area_Ave )

end

fprintf('\n\n\n')

end



