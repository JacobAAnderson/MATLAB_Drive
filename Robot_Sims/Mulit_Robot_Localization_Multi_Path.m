% Multi-robot localization Sim -- Stats
% Jacob Anderson
% 1/30/2020

close all
clear all
clc
format long g

ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(false);                                                     % Indicate whether new data should be selescted by the user


numPaths = 2; % 50
trials = 1;   % 100
num_auv = 2;


%% Get Data Sources
file_bathy = ui.GetFile('bathymetry','mat');                                % Get file path GUI
if isempty(file_bathy), return, end                                         % End script if there isn't an input for it to use

for jj = num_auv: -1: 1
    file_wp{jj} = ui.GetFile( sprintf('waypoints%d', jj),'txt');             % Get file path GUI
end

if any(isempty(file_wp)), return, end

file_saveResults = ui.SaveFile('Results','mat');                            % Get File Path to save Results
if isempty(file_saveResults), return, end                                   % End script if there isn't an input for it to use

clear ui jj;


%% Computer parimiter
bathy = LoadBathymetry(file_bathy);                                         % Load Bathymetry Map
bathy = bathy.ReduceResolution(0.5);

perim = GetPerim(bathy);

% fig = bathy.Plot_3DModel(0, 90);
% hold on
% plot(perim(:,1),perim(:,2), '*')



%% Track Simulation Results
count = 1;
errorCount = 1;

simulations = {'DR', 'PF_EKF', 'PF_TwoWay_Acoms', 'PF_ConstraintOnly_Acoms', 'Policy', 'rand05', 'rand1', 'rand2'};
dataTypes   = {'Err1', 'Err2', 'Time' 'Comms1', 'Coms2'};


rAd_master = RandomAssData([3000, trials*numPaths], simulations, dataTypes);    % Instanciate RandomAssData

for p = 1 : numPaths, fprintf( '\n\nPath number %d of %d\n\n', p, numPaths)
    
%     try
        for vv = num_auv:-1:1
            path = MakePath(perim);
%             writematrix(path(vv,:), file_wp{vv})
            dlmwrite(file_wp{vv}, path, ',') 
        end
        
%         rAd_path = Simulation1(trials, simulations, bathy, file_saveResults, file_wp1, file_wp2, count, false);
        rAd_path = Simulation2(trials, num_auv, simulations, bathy, file_saveResults, file_wp, count, true);
        
        rAd_master = rAd_master.AddRad(rAd_path);
        
        save('20200304_rAd_Master.mat', 'rAd_master')
        
        count = count+1;
        
%     catch
%         errorCount = errorCount + 1; 
%     end
    
end

rAd = rAd.Eval_Stats;
save('20200304_rAd_Master.mat', 'rAd_master')


% Show Individual AUV Results
figure('name', 'Individual Error')
subplot(2,1,1)
rAd_master.PlotStat(simulations, 'Ave', 'Err1', 'Time')

subplot(2,1,2)
rAd_master.PlotStat(simulations, 'Ave', 'Err2', 'Time')


% Show Joint AUV Results
figure('name', 'Joint Error')
rAd_master.PlotStat(simulations, 'Ave', 'JointErr', 'Time')




fprintf('\n\n\nSimulation errored out %d times\n\n',errorCount)



%% Functions
function perim = GetPerim(bathy)

map = bathy.AsMapp;
perim = map.Perim;

for ii = 1:5
    ind = 2:2:size(perim,1);
    perim(ind,:) = [];
end

perim = simplifyPolygon(perim, 1);
loops = expandPolygon(perim, -200);

sz = cellfun(@(x) size(x,1), loops);

perim = loops{sz == max(sz)};

end


function path = MakePath(perim)

tf = true;
while tf
    
    X = randi(size(perim,1), 2);
    
    path = [perim(X(1),:), 0; perim(X(2),:), 0];
    
    V = path(1,1:2) - path(2,1:2);
    D = vecnorm(V);
    V = V/D;
    
    % Check for pathes leaving the map and / or in collitions
    ii = 0 : 10: D;
    test = path(2,1:2) + ii'.*V;
        
    if all( inpolygon(test(:,1), test(:,2), perim(:,1), perim(:,2)) ) && D > 500
        tf = false;
    end
    
end
end






