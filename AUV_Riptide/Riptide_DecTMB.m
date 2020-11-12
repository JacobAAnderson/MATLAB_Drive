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

timeLine_mission = [];                                                      % Keep track of which time line we have


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

disp('Saving Riptide Data')
save("/Users/jake/_OutPuts/RT_FostersLake.mat", 'RT', 'geotiff', 'bathy')

clear owtof filters



%% Select Mission
close all
clc

load('/Users/jake/_OutPuts/RT_FostersLake.mat')

% ------ 2020-09-21 ----------------
% missions = [ 1, 1];         % Fosters Lake Mission Set: CurcumNav1, CircumNav2
% missions = [ 2, 2];         % Fosters Lake Mission Set: CurcumNav1, CircumNav2        --> Good Acoms
% missions = [ 3, 3];         % Fosters Lake Mission Set: Dimond1_A, Dimond1_B          --> Decent Acomms

% ------ 2020-10-09 ----------------
% missions = [ 5, 7];         % Fosters Lake Mission Set: test1, test2
% missions = [ 7, 8];         % Fosters Lake Mission Set: ZZ1, ZZ2
% missions = [ 8, 9];         % Fosters Lake Mission Set: ZZ1, ZZ2
% missions = [ 9,10];         % Fosters Lake Mission Set: ZZ3_A, ZZ3_B                  --> Decent Acoms
% missions = [10,11];         % Fosters Lake Mission Set: ZZ3_A, ZZ3_B                  --> Decent Acomms

% ------ 2020-10-27 ----------------
 missions = [11,12];         % Fosters Lake Mission Set:CircumNav2_B, CircumNave2_A    --> Decent Acomms
% missions = [14,15];         % Fosters Lake Mission Set:ZZ1, ZZ2                         --> Decent Acomms
% missions = [15,16];         % Fosters Lake Mission Set:Dimond1_A, Dimond1_B  (DEC-TBN Does better)
% missions = [16,17];         % Fosters Lake Mission Set:ZZ4_B, ZZ4_A  



% fig1 = geotiff.Show;                                                      % Show geotiff for plotting paths

for v = 1: numel(RT)
    
    % -- Select mission --
    RT(v) = RT(v).Select_Mission( missions(v) );                            % Choose Mission
    RT(v) = RT(v).Get_VehiclePaths;                                         % Evaluate GT and DR paths for speed and compass modeling
    
    % -- Terrain Based Navigation stuff --
    RT(v) = RT(v).Add_TBN(bathy);                                           % Prep Particle filters
    RT(v) = RT(v).Model_VehicleSpeed;                                       % Calculate corrected speed from GPS
    RT(v) = RT(v).Model_Compass('poly2');                                   % Model / Calibrate Compass data
    RT(v) = RT(v).Model_Altimiter(false);                                   % Filter Altimiter data and build normal distribution
    RT(v) = RT(v).Get_VehiclePaths;                                         % Evaluate GT and Corrected DR paths from the speed and compass models
    
    
    % -- Plot/display Stuff --
    RT(v).Disp_Manifest;                                                    % Display Vehicle Manifest
%     fig1 = RT(v).Plot_Paths(fig1);                                        % Show paths on a geotiff
%     RT(v).Plot_AltimeterProfile;                                          % Show Altimeter profile
%     RT(v).Plot_Altimeter(geotiff);                                        % Show Altimeter location on map

end





[mu, sig] = rtAcms.Mission_Latency(RT);                                     % Get Acoustic modem model for the selected mission
                                   

for v = 1:numel(RT)                                                         % Set Acomms Noise Distributions 
    m = seconds(mu(v));
    s = seconds(sig(v))*sos;
    RT(v).tbn = RT(v).tbn.SetNoise('acoms',  'Normal', m, s * 1.5);         % Set acoustic modem error Distributions
end

rtAcms.Plot_AcousticCommunications(RT, geotiff);                            % Plot Acoustic Communications

clear v m s

RT_ = RT;                                                                   % Keep a clean copy of RT


% ---- Get Timeline -----
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



%% Communication planning
clc
close all

for ii = [1: numel(RT); numel(RT):-1:1], v = ii(1); s = ii(2);

    RT(v) = RT(v).Add_CommsPlanner;                                                     % Instanciate communications planner
    
    RT(v).commsPlanner = RT(v).commsPlanner.Add_Vehicle(RT(v).VehicleInfo('self'));     % Load vehicle info into the planner
    RT(v).commsPlanner = RT(v).commsPlanner.Add_Vehicle(RT(s).VehicleInfo);
   
end

coms_times = RT(1).acomms.gpsDate;

dt = seconds(mean ( coms_times(2:end) - coms_times(1:end-1) ));

sz = numel(coms_times);
coms_policy = zeros(numel(RT), sz);

for v = 1: numel(RT)
    fprintf('Comms Planning for %s\n', RT(v).name)
    
    policy = RT(v).commsPlanner.Plan1(sz, dt, 1.25, true);
    
    coms_policy(v, 1:size(policy,1) ) = policy';
end


clear ii v s


%% _____ Do TBN _____________________
clc
close all

num_trials = 1;

% --- Instanciate rAd to track data ---
paths = {'DR_Correctd', 'DEC_TBN', 'TBN'};                                  % Set up rAd with DR corrected data
types = {'Dory_Error', 'Nemo_Error', 'Joint_Error', 'Time'};

[rAd, gtPaths] = GetRAD(RT, paths, types, num_trials, timeLine);            % Create rAd with DR Correctd data entered

RT = RT_;


for prop = { 'altimeter', 'compass' 'speed';           % --> Type of noise
             'Normal',    'Normal', 'Normal';          % --> Default distribion type
              0,           0,        0;                % --> Default mu
              2.5,         30,       2.0;              % --> Sigma for Dory
              2.5,         30,       2.0}              % --> Sigma for Nemo
    

        RT(1).tbn = RT(1).tbn.SetNoise(prop{1}, prop{2}, prop{3}, prop{4});
        RT(2).tbn = RT(2).tbn.SetNoise(prop{1}, prop{2}, prop{3}, prop{5});

end



% --- TBN Simulation ---
for dec_tbn = [false, true]                                                 % Switch between TBN and DEC-TBN
    for iter = 1: num_trials
        
        [tbnPaths, sim_times] = Riptide_TBN_Sim(RT, timeLine, sos, mu, dec_tbn);    % Do TBN / DEC-TBN simulation
        
        rAd = TallyPathData(rAd, iter, paths(2:3), gtPaths, tbnPaths, sim_times, timeLine, dec_tbn);     % Calculate path errors and add to rAd
        
        
        if iter == num_trials
            
            for v = 1: numel(RT)
                
                if dec_tbn
                    RT(v) = RT(v).Add_Path(tbnPaths(:,:,v), "DEC_TBN", 'utm');
                else        
                    RT(v) = RT(v).Add_Path(tbnPaths(:,:,v), "TBN",     'utm');
                end
                
            end
        end
          
    end
end


rAd = rAd.Eval_Stats('Time');                                               % Determine average error and total error

DisplyRAD(rAd, types, paths)                                                % Display Result Stats


% ____ Display Specific Graphs ___________________________________________
for v = 1:numel(RT)
    
    RT(v).Plot_Paths(geotiff, {'GT', 'DR_corrected', 'TBN', 'DEC_TBN'});                                             % Plot lat-lon path on geotiff
   
end





%% ========= Functions ==============


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
    drPaths(1:sz,:,v) = RT(v).path.DR_corrected.utm;
    
    drTimes(1:sz,v) = RT(v).data.gpsDate;
    
    
end

rAd = TallyPathData(rAd, 1, paths(1), gtPaths, drPaths, drTimes, timeLine, true);


end




%% Add Path data to rAd 
function rAd = TallyPathData(rAd, iter, paths, gtPaths, tbnPaths, times, t, dec_tbn)

if dec_tbn, path = paths(1);
else,       path = paths(2);
end

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
    rAd = rAd.Add_Result(iter, path{1}, err(ii,1), err(ii,2), err(ii,3), t(ii) ); % Add data to the respective set
end


rAd = rAd.ResetIndex;


end



%% Display rAd stats
function DisplyRAD(rAd, types, paths)

%for t = types(1:3), rAd.PlotData( paths, t{1}, 'Time'), end
for t = types(1:3), rAd.PlotStat( paths, 'Ave', t{1}, 'Time'), end


fprintf('\n\n\n')
fprintf("Total Joint Error for DR Cor.:  %.2f [m s]\n", rAd.Stats.DR_Correctd.Joint_Error.Area_Ave)
fprintf("Total Joint Error for DEC-TBN:  %.2f [m s]\n", rAd.Stats.DEC_TBN.Joint_Error.Area_Ave)
fprintf("Total Joint Error for TBN:      %.2f [m s]\n", rAd.Stats.TBN.Joint_Error.Area_Ave)
fprintf('\n\n\n')

end



