% Demo Bathymetry_Map
% Jacob Anderson
% 9/16/2019

close all
clear all
clc


%% Choose Data Type --> Uncomment the block that you want

%___ Load Ecomapper Mission Logs ________________________
fileType = 'mat';
fun = @(file) getData(file);
% fun = @(file) Logdoc2mat(file);


%% Load Data
ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(true);                                                     % Indicate whether new data should be selescted by the user

filePath = ui.GetDir('Data',fileType);                                      % Get folder GUI

if isempty(filePath)                                                        % Check if the Input box was cnaceled
    disp('Get Directory Input box Cancelled')
    return                                                                  % End script if there isn't an input for it to use
end

Data = BatchDir(filePath,fileType,fun);


if isempty(Data)                                                            % Check that data is present
    disp("Empty Data Structure")
    disp("Doulbe check the directroy that you chose")
    return                                                                  % End script if there isn't any data to process
end

A = cat(1,Data.LNH_regionData);

%% Make Map from other data
bathy_map = Bathymetry_Map;
bathy_map = bathy_map.MakeMap(A);

% Access Bathymetry Data
depth    = bathy_map.Elevation;                                             % Elevation Profile
sigma    = bathy_map.Variance;                                              % Uncertainty in the elevation
XX       = bathy_map.Easting;                                               % X utms
YY       = bathy_map.Northing;                                              % Y utms
utmZone  = bathy_map.UTM_Zone;                                              % Utm Zone
dataLogs = bathy_map.Logs;                                                  % Data lgs used to create the map
mapInfo  = bathy_map.Header;                                                % Discription of the data fields


% Display Bathymetry
fig1 = bathy_map.Plot_3DModel(250, 30);                                     % Plot Bathymetry Model



% Save Map
saveFile = ui.SaveFile('Bathymetry','mat');                                 % GUI to get file path to save variable / object
if isempty(saveFile), return, end                                           % End script if there isn't an input for it to use

save(saveFile, 'bathy_map')


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


function data = getData(file)

    data = load(file);

end

