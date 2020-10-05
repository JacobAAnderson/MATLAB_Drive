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
ui = ui.NewData(false);                                                      % Indicate whether new data should be selescted by the user


% Get file with file path GUI ----------------------------------------------------------------------------------------------
path_data = ui.GetDir('Data','csv');                                        % GUI to get path to data logs
if isempty(path_data), return, end                                          % End script if there isn't an input for it to use

path_acomms = ui.GetDir('Acomms','csv');                                    % GUI to get file path to acoustic communication logs
if isempty(path_acomms), return, end 

path_pAcommsHand = ui.GetDir('pAcomms','txt');                              % GUI to get file path to pAcommsHandler logs
if isempty(path_pAcommsHand), return, end 

path_waypoints = ui.GetDir('waypoints','txt');                              % GUI to get file path mission waypoint courses
if isempty(path_waypoints), return, end

geotiff = ui.GetGeoTiff;                                                    % GUI to get Geotiff
if isempty(geotiff.Image), return, end 


%% Instanciate Riptide Data and Add data
filters = {'GPS_fix', 0};                                                   % Data filter {"Data_Field", value to be filterred out}
%           'ALT_ALTITUDE', 0};        

RT = Riptide_Data(path_data, '210', filters);                               % Create Object with Riptide Data data logs and give vehicle name '210'

RT = RT.Add_Acomms(path_acomms, filters);                                   % Add acoustic communications logs
RT = RT.Add_pAcommsHandler(path_pAcommsHand);                               % Add pAcommsHandler logs
RT = RT.Add_Waypoints(path_waypoints);                                      % Add mission waypoints

RT = RT.Get_Manifest("all");                                                % Create a table summerising the data and acomms logs, skip instances where the vehicle is idel and disply the table

clear path_data path_acomms path_pAcommsHand path_waypoints ui              % Clean up the work space


%% Do Stuff with data
close all
clc
RT.Disp_Manifest;                                                           % Print manifest

RT = RT.Select_Mission(5);                                                  % Select a mission

RT = RT.Model_VehicleSpeed;                                                 % Determine the true speed of the vehicle
RT = RT.Get_VehiclePaths;                                                   % Get GPS and dead reckoning paths

fig1 = RT.Plot_Paths(geotiff);                                              % Show paths on a geotiff
fig1 = RT.Plot_Course(fig1);                                                % Show mission waypoints on the figure
fig1 = RT.Plot_Acomms(fig1);                                                % Show Acoustic communications



