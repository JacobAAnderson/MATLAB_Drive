% Demo Side Scan Sonar Data Class
% Jake Anderson
% 9/26/2019

close all
clear all
clc


ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(false);                                                     % Indicate whether new data should be selescted by the user

filePath = ui.GetDir('SideScanSonarData','logdoc');                         % Get Side Scan Sonar .logdoc file
if isempty(filePath), return, end                                           % End script if there isn't an input for it to use


sssData = SSS_Data(filePath);