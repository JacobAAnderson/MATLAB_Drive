% Demo Riptide_AcommsLog2Mat
% Jacob Anderson
% Robotic Decition Making Laboratory (RDML)
% Oregon State University
% Corallis OR
% August 20, 2020

close all
clear all
clc


ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(true);                                                     % Indicate whether new data should be selescted by the user


% Get file with file path GUI ----------------------------------------------------------------------------------------------
file1 = ui.GetFile('Acomms1','csv');                                        % Get file path GUI
if isempty(file1), return, end                                              % End script if there isn't an input for it to use

% file2 = ui.GetFile('Acomms2','csv');                                        % Get file path GUI
% if isempty(file2), return, end                                              % End script if there isn't an input for it to use

ui = ui.NewData(false);                                                     % Indicate whether new data should be selescted by the user
geotiff = ui.GetGeoTiff;
if isempty(geotiff.Image), return, end                                      % End script if there isn't an input for it to use


Data1 = Riptide_AcommsLog2Mat(file1);
% Data2 = Riptide_AcommsLog2Mat(file1);

clear file1 file2 ui

%% Do stuff with the data
close all


lat1  = Data1.vehicle(:,1);
lon1  = Data1.vehicle(:,2);

figure('Name',"Map", 'numbertitle','off')
geotiff.Show
hold on
plot(lon1, lat1, '*-w')
% plot(lon2, lat2, '*-b')


% --- Show Acoustic Communications ---
% [~] = Plot_AcousticCommunications(RT, geotiff);




