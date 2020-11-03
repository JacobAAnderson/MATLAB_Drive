% Demo Riptide_Data
% Jacob Anderson
% Robotic Decition Making Laboratory (RDML)
% Oregon State University
% Corallis OR
% August 25, 2020

close all
clear all
clc

ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(false);                                                     % Indicate whether new data should be selescted by the user

% name = 'Dory';
name = 'Nemo';

% Get file with file path GUI ----------------------------------------------------------------------------------------------
paths.data = ui.GetDir( sprintf('%s Data',name), 'csv');                    % GUI to get path to data logs
if isempty(paths.data), return, end                                         % End script if there isn't an input for it to use

paths.acomms = ui.GetDir( sprintf('%s Acomms',name), 'csv');                % GUI to get file path to acoustic communication logs
if isempty(paths.acomms), return, end 

paths.pAcommsHand = ui.GetDir( sprintf('%s pAcomms', name), 'txt');         % GUI to get file path to pAcommsHandler logs
if isempty(paths.pAcommsHand), return, end 


paths.waypoints = ui.GetDir('Waypoints','txt');                             % GUI to get file path mission waypoint courses
if isempty(paths.waypoints), return, end

paths.offset = ui.GetFile('Off Set','txt');                                 % GUI to get Get water level offset

geotiff = ui.GetGeoTiff;                                                    % GUI to get Geotiff
if isempty(geotiff.Image), return, end 


%% Instanciate Riptide Data and Add data
filters = {'GPS_fix', 0};                                                   % Data filter {"Data_Field", value to be filterred out}
%           'ALT_ALTITUDE', 0};        

RT = Riptide_Data(paths.data, paths.offset, name, filters);               % Create Object with Riptide Data data logs and give vehicle name '210'

RT = RT.Add_Acomms(paths.acomms, filters);                                  % Add acoustic communications logs
%RT = RT.Add_pAcommsHandler(paths.pAcommsHand);                              % Add pAcommsHandler logs
%RT = RT.Add_Waypoints(paths.waypoints);                                     % Add mission waypoints

RT = RT.Get_Manifest("all");                                                % Create a table summerising the data and acomms logs, skip instances where the vehicle is idel and disply the table

clear paths ui                                                              % Clean up the work space


%% Select Mission
close all
clc
RT.Disp_Manifest;                                                           % Print manifest

RT = RT.Select_Mission(17);                                                 % Select a mission
RT = RT.Model_VehicleSpeed;                                                 % Determine the true speed of the vehicle
RT = RT.Model_Compass('poly2');                                             % Model / Calibrate Compass data
RT = RT.Get_VehiclePaths;                                                   % Get GPS and dead reckoning paths

fig1 = RT.Plot_Paths(geotiff);                                              % Show paths on a geotiff
% fig1 = RT.Plot_Course(fig1);                                                % Show mission waypoints on the figure
% fig1 = RT.Plot_Acomms(fig1);                                                % Show Acoustic communications



%% Model / Calibrate Altimiter Data
fig = RT.Plot_Altitude;                                                     % Show Altimeter profile

RT = RT.Model_Altimiter(20);                                                % Filter Altimiter data and build normal distribution
[~] = RT.Plot_Altitude(fig);                                                % Show Altimeter profile

% [~] = RT.Plot_Alt_Geotiff(geotiff);                                         % Show Altimeter location on map
    



