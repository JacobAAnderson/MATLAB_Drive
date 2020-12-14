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


%% Access Bathymetry Data
% depth    = bathy_map.Elevation;                                             % Elevation Profile
% sigma    = bathy_map.Variance;                                              % Uncertainty in the elevation
% XX       = bathy_map.Easting;                                               % X utms
% YY       = bathy_map.Northing;                                              % Y utms
% utmZone  = bathy_map.UTM_Zone;                                              % Utm Zone
% dataLogs = bathy_map.Logs;                                                  % Data lgs used to create the map
% mapInfo  = bathy_map.Header;                                                % Discription of the data fields


%% Display Bathymetry
fig1 = bathy_map.Plot_3DModel(250, 30);                                     % Plot Bathymetry Model with view angle

if ~isempty(geotiff.Image)
    fig2 = bathy_map.PlotMap(geotiff);                                      % Plot Map over the Geotiff 
end


%% Save Map
ui = ui.NewData(true);                                                      % Indicate that new data should be selescted by the user
saveFile = ui.SaveFile('Bathymetry','mat');                                 % GUI to get file path to save variable / object
if isempty(saveFile), return, end                                           % End script if there isn't an input for it to use
disp("Saving Bathymetry Map")
save(saveFile, 'bathy_map')


%% Other Properties
bathy_map = bathy_map.ReduceResolution(0.5);                                % Reduce Bathymetry map to 50% Resolution

mapp = bathy_map.AsMapp;                                                    % Get Bathymetry as a Mapp


%% Export Bathymetry Map
% saveFile = ui.SaveFile('Bathymetry','stl');                                 % GUI to get file path to save variable / object
% if isempty(saveFile), return, end                                           % End script if there isn't an input for it to use
% 
% bathy_map.Eport(saveFile);
% 
% % Examin Bathymetry STL
% [x, y, z, c] = stlread(saveFile);
% 
% figure('NumberTitle','off','Name','Origonal Object')                        % Show what came in
% axis equal
% patch(x, y, z, c, 'FaceAlpha', 1)
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% view(45,30)
% 


