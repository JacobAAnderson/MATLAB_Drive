% Riptide AcommsHandler output parser
% Jacob Anderson
% Robotic Decition Making Laboratory (RDML)
% August 28, 2020

close all
clear all
clc

format long g

ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(true);                                                      % Indicate whether new data should be selescted by the user

% Get file with file path GUI ----------------------------------------------------------------------------------------------
file = ui.GetFile('Data','txt');                                            % Get file path GUI
if isempty(file), return, end                                               % End script if there isn't an input for it to use

acomms = Riptide_pAcommsHandler2Mat(file);