% Multi-robot localization Sim -- Stats
% Jacob Anderson
% 1/30/2020

close all
clear all
clc
format long g

ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(false);                                                     % Indicate whether new data should be selescted by the user

num_auv = 2;
trials = 1;
simulations = {'PolicyNo', 'PolicyYes'};
maxErr = 1.25;

%% Get Data Sources
file_bathy = ui.GetFile('bathymetry','mat');                                % Get file path GUI
if isempty(file_bathy), return, end                                         % End script if there isn't an input for it to use

for jj = num_auv: -1: 1
    file_wp{jj} = ui.GetFile( sprintf('waypoints%d', jj),'txt');            % Get file path GUI
end

if any(isempty(file_wp)), return, end

file_saveResults = ui.SaveFile('Results','mat');                            % Get File Path to save Results
if isempty(file_saveResults), return, end                                   % End script if there isn't an input for it to use

clear ui


%% Track Simulation Results

%--> Simulations: 'PF_EKF' 'PF_OneWay_Acoms' 'PF_TwoWay_Acoms'  'PF_OneWay_Acoms_GT' 'PF_ConstraintOnly_Acoms'  'Policy' 'rand05'  'rand1'  'rand2'

bathy = LoadBathymetry(file_bathy);                                         % Load Bathymetry Map
bathy = bathy.ReduceResolution(0.5);

% rAd = Simulation1(trials, simulations, bathy, file_saveResults, file_wp1, file_wp2, false);
% rAd = Simulation2(trials, num_auv, simulations, bathy, file_saveResults, file_wp, 1, false);
[rAd, coms_policy] = Simulation3(trials, num_auv, maxErr, simulations, bathy, file_saveResults, file_wp, 1, true);


clear bathy trials


%% Show Individual AUV Results

% figure('name', 'Individual Error')
% subplot(2,1,1)
% rAd.PlotStat(simulations, 'Ave', 'Err1', 'Time')
% 
% subplot(2,1,2)
% rAd.PlotStat(simulations, 'Ave', 'Err2', 'Time')


%%  Show Joint AUV Results


% simulations = {'Policy', 'PF_EKF', 'half', 'quart'};

rAd.ErrorPlotStat(simulations, 'Ave', 'JointErr', 'Time', 'SEM')



