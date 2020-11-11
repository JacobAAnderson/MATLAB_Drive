
% Demo RandomaAssData
% Jake Anderson
% 12/6/2019

close all
clear all
clc


%% File Pathes
ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(false);                                                     % Indicate whether new data should be selescted by the user

saveResults = ui.SaveFile('Results','mat');                                 % Get File Path to save Results
if isempty(saveResults), return, end                                        % End script if there isn't an input for it to use

saveStats = ui.SaveFile('Stats','tex');                                     % Get File Path to save statistics from results
if isempty(saveStats), return, end                                          % End script if there isn't an input for it to use


%% Create Data
trials = 20;                                                                % Number of trial runs
rAd = RandomAssData(trials, {'env1'; 'env2'; 'env3'}, {'Err'; 'Perc'; 'Time'});    % Instanciate RandomAssData

for iter = 1 : trials
    
    fprintf("Iteration: %i\n", iter)
    
    for env = 1:3                                                           % Cycle through the data environments
        
        mse      = rand;                                                    % Data to be added to the data set
        per      = rand;            
        sim_time = rand;
        
        switch env
            case 1, rAd = rAd.Add_Result('env1', mse, per, sim_time);       % Add data to the respective set
            case 2, rAd = rAd.Add_Result('env2', mse, per, sim_time);
            case 3, rAd = rAd.Add_Result('env3', mse, per, sim_time);
            otherwise
                continue
        end
        
    end
end

save(saveResults, 'rAd')                                                    % Save Data


%% Do statisics on the data and output them as a table
rAd = rAd.Eval_Stats;                                                       % Do statistics 

rAd.Write_Stats(saveStats, 'Erro', '%1.5f')                                 % Write statistics to a text file as a table
rAd.Write_Stats(saveStats, 'Perc', '%1.4f')

