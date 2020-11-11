
close all
clear all
clc

% load('/home/jake/MATLAB-Drive/20200304_rAd_Master.mat')
load('/Users/jake/Downloads/20200304_rAd_Master.mat')

% flds = fields(rAd.Data);
% 
% for ii = 1: numel(flds) 
%     err = rAd.Data.(flds{ii}).Err1 + rAd.Data.(flds{ii}).Err2;
%     rAd.Data.(flds{ii}).JointErr = err;
% end

%%
rAd_master = rAd_master.Eval_Stats;

%%
% simulations = fields(rAd.Data);
% simulations = { 'DR', 'PF_EKF', 'PF_OneWay_Acoms', 'PF_TwoWay_Acoms', 'rand05', 'rand1', 'rand2', 'Policy'};

% simulations = {'DR', 'PF_EKF', 'PF_TwoWay_Acoms', 'rand1', 'Policy'};
simulations = {'PF_EKF', 'rand05', 'rand1', 'rand2', 'Policy'};

% Show Individual AUV Results
figure
subplot(2,1,1)
rAd_master.PlotStat(simulations, 'Ave', 'Err1', 'Time')

subplot(2,1,2)
rAd_master.PlotStat(simulations, 'Ave', 'Err2', 'Time')

% Show Joint AUV Results
figure
rAd_master.ErrorPlotStat(simulations, 'Ave', 'JointErr', 'Time')

