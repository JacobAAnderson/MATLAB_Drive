% Demo Script for UserInputs
% Jacob Anderson
% 1/1/2020

close all
clear all
clc


ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(true);                                                     % Indicate whether new data should be selescted by the user


% % Get folder GUI ------------------------------------------------------------------------------------------------------
% filePath = ui.GetDir('Data','pdf');                                         % Get folder GUI
% if isempty(filePath), return, end                                           % End script if there isn't an input for it to use
% 
% 
% % Get file with file path GUI ----------------------------------------------------------------------------------------------
% file = ui.GetFile('Data','pdf');                                            % Get file path GUI
% if isempty(file), return, end                                               % End script if there isn't an input for it to use
% 
% 
% % Save file GUI ---------------------------------------------------------------------------------------------------------
% saveFile = ui.SaveFile('info','txt');                                       % GUI to get file path to save variable / object
% if isempty(saveFile), return, end                                           % End script if there isn't an input for it to use
% 
% 
% % Get Geotiff GUI ------------------------------------------------------------------------------------------------------
% lon = [-107.9113366574420,  -107.9073118934940];                            % Min and Max longitudes
% lat = [  37.2387203649409,    37.2387286862432];                            % Min and Max latitudes

% 44°29'08.9"N 122°25'52.3"W

geotiff = ui.GetGeoTiff;                                                    % Get Geotiff
% geotiff = ui.GetGeoTiff( lon, lat );                                      % Gui to get Geotiff ---> (lon, lat) is optional for locating file on your computer but needed for retriving map from Google / web map server
if isempty(geotiff.Image), return, end                                      % End script if there isn't an input for it to use

geotiff.Show


