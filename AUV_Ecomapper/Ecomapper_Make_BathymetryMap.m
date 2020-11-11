% Demo Bathymetry_Map
% Jacob Anderson
% 9/16/2019

close all
clear all
clc


% Get Data Sources
ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(false);                                                     % Indicate whether new data should be selescted by the user

geotiff = ui.GetGeoTiff;                                                    % Gui to get Geotiff ---> (lon, lat) is optional for locating file on your computer but needed for retriving map from Google / web map server


%% Bathymetry Map from Ecommaper data
filePath = ui.GetDir('Ecomapper','log');                                    % Get folder GUI
if isempty(filePath), return, end                                           % End script if there isn't an input for it to use

geotiff = ui.GetGeoTiff;                                                    % Gui to get Geotiff ---> (lon, lat) is optional for locating file on your computer but needed for retriving map from Google / web map server

% Make Bathymetry map from Ecomapper logs
station_ID = 9410079;                                                       % NOAA station ID for Catlina IS tide data
bathy_map = Bathymetry_Map(filePath, station_ID);                           % Create Bathymetry Map uisng NOAA tide data


%% Display Bathymetry
fig1 = bathy_map.Plot_3DModel(250, 30);                                     % Plot Bathymetry Model

if ~isempty(geotiff.Image)
    fig2 = bathy_map.PlotMap(geotiff);                                      % Plot Map over the Geotiff 
end


%% Save Map
ui = ui.NewData(true);                                                      % Indicate that new data should be selescted by the user
saveFile = ui.SaveFile('Bathymetry','mat');                                 % GUI to get file path to save variable / object
if isempty(saveFile), return, end                                           % End script if there isn't an input for it to use
disp("Saving Bathymetry Map")
save(saveFile, 'bathy_map')


