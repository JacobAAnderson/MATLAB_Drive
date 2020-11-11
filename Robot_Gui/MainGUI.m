%% Robotics GUI.m
% Jacob Anderson
% 12-01-2017

% This is a program that generates a GUI for doing environmental niche modeling robot navigation and More.
% 12/2017 --> GUI window class created
% 01/2018 --> Environmental Nich Modeling tool created
% 02/2018 --> Ocupancy Grid and Georeferancing tool created
% 07/2019 --> Side Scan Sonar tool created
% 09/2019 --> Side Scan Feature Predition Analysis tool Created

clc, close all, clear all

global ui               % Variable to save and access saved user data
% global DATA             % Data from vehicle
% global Model            % Area model being developed by this app

% global MainWindow       % Main GUI window
% global SSSWindow        % Side Scan Soanr Window
% global EnvWindow        % Environmental Niche Modeling GUI window

[filepath,name,~] = fileparts(mfilename('fullpath'));   % Find path to this file

ui = SavedUserInputs(name);                             % Create an file for saving data or acces previously saved data
ui = ui.NewData(true);                                  % Indicate whether new data should be selescted by the user

% Make sure subfolder are included in search path
toolbox = fullfile(filepath,'toolbox');
source  = fullfile(filepath,'src');
addpath(toolbox,source);

% Get a list of tools in the toolbox
toolbox = dir(toolbox);                                 % Get list if file in the toolbox directory
toolbox = toolbox(~cat(1,toolbox.isdir));               % Filter out invalid entries
toolbox = {toolbox.name};                               % Isolate file manes

% Remove file extention and check that the file extention are appropriate (.m)
jj = 1;
tools(length(toolbox)) = {'Tools in tool box'};         % Prealocate memory 

for ii = 1:length(toolbox)                              % Cycle through the files in the tool box and seperate file extention
    str = strsplit(toolbox{ii}, '.');
    
    if contains(str(2),'m~'), continue                  % These files are backup for files open in the matlab editor, don't add them to the tool box
    elseif contains(str(2),'m')                         % Check that the file extention is right
        tools(jj) = str(1);
        jj = jj + 1;
    end
end

tools(jj:end) = [];                                     % Clean up extra entries there were prealocated for memeory

Start(tools)                                            % Start the GUI

clear filepath name toolbox tools source jj ii str

%% GUI Window
function Start(tools)
global MainWindow   % Main GUI window

windowName = 'Marine Robotic Dev GUI';
panelName = 'Main';

MainWindow = GUIwindow(windowName);                                                             % Instanciate the GUI
MainWindow.GUIpanels(panelName,'make','position',[0.8 0 0.2 1],'backgroundcolor',[0, 0.4470, 0.7410]);  % Create putton panel for primary GUI functions

MainWindow.Header.String = 'Select Data Input';                                                 % Prompt to select data

% Make buttons in the button panel --> { Buttin type,  Button Name,  Callback function and their inputs,  Callback function for slider listeners }
MainWindow.GUIbuttons( panelName,'add', ...
   {'pushbutton','Load Geotiff',      {@GetGeoTiff,     panelName},[], 'Select a goe-referances image\nFile extensions: .tif, .tiff';                                                                                       %  1
    'pushbutton', 'Color Threshold',  @ColorThreshold,             [], 'Color threshold Geotiff\nColor threshold is used to seperate land from water\nUse "Make Occupancy Grid function to apply the color thresholding';   %  2
    'space',      ' ',                [],                          [], '';                                                                                                                                                  %  3
    'pushbutton', 'Load Shape File',  {@GetShapeFile,   panelName},[], 'Select a list of waypoints that define the area of interest\nFile extensions: .txt';                                                                %  4
    'pushbutton', 'Clear Shape File', {@ClearShapeFile, panelName},[], '';                                                                                                                                                  %  5
    'space',      ' ',                [],                          [], '';                                                                                                                                                  %  6
    'pushbutton', 'Load Data',        {@LoadData,       panelName},[], 'Load Sensor Data\nChose a directory containing the data files\n All files in that directory will be read in\nFile extensions: .log';                %  7
    'pushbutton', 'Clear Data',       {@ClearData,      panelName},[], '';                                                                                                                                                  %  8
    'space',      ' ',                [],                          [], '';                                                                                                                                                  %  9
    'text'        'Tool Box',         [],                          [], '';                                                                                                                                                  % 10
    'popupmenu', tools,               {@OpenToolbox,    panelName},[], '';                                                                                                                                                  % 11
    'space',      ' ',                [],                          [], '';                                                                                                                                                  % 12
    'pushbutton', 'Clear All',        @ClearAll,                   [], '';                                                                                                                                                  % 13
    'pushbutton', 'Save',             @SaveAll,                    [], '';                                                                                                                                                  % 14
    });


% Hide buttons that are not ready of use
Panel = MainWindow.Panels(panelName);   % Put the new version of the panel into its mapped container
Panel.button(2).Push.Visible = 'off';  % Color Threashold
Panel.button(5).Push.Visible = 'off';  % Clear shape file
Panel.button(8).Push.Visible = 'off';   % Clear Data 

end


%% Load Data
function LoadData(~,~,panelName)
global MainWindow ui Model DATA

MainWindow.Header.String = 'Select Folder Containing Data Files';           % Prompt user on what to do

filePath = ui.getDir('EM_Data','log');                                      % Get folder GUI
if isempty(filePath), return, end                                           % End script if there isn't an input for it to use

MainWindow.Header.String = 'Loading Data, Please Wait';                     % Prompt user on what to do
set(MainWindow.Figure, 'pointer', 'watch')
drawnow;

DATA.EM_data = EM_Data(filePath);                                           % Run Function to import Ecomapper data

% % Get file and folder GUI ----------------------------------------------------------------------------------------------
% [fileName, filePath] = ui.getFile('WaterOffSet','txt');             % GUI to get the name and file path of a file
% fprintf('\nFile Path and name: %s %s \n\n', filePath, fileName);
%
% if fileName == 0                                                            % Check if the Data Input box was cnaceled
%     disp('Input box Cancelled')
%     waterOffSetLog = '';
% else
%     waterOffSetLog = fullfile(filePath,fileName);                           % Get full file path to water offset log
% end

if isfield(Model,'Boundaries')
    
    answer = questdlg('Apply Boundary to the Data?', ...
        'Data Filter', ...
        'Yes','No','No');

    switch answer
        case 'Yes'
            filter = @inpolygon;
            DATA.EM_Data = applyfilter(DATA.EM_Data, 'vehicle', filter, Model.Boundaries);        
    end
end

dataPoints = cat(1,DATA.EM_data.RawData.vehicle);

set(MainWindow.Figure, 'pointer', 'arrow')
drawnow;

MainWindow.Addlayer('DataPoints', dataPoints, panelName);                       % Add datapoints to the GUI

MainWindow.Header.String = 'Make Ocupancy Grid or Environmantal Niche Model';   % Prompt to select data

Panel = MainWindow.Panels(panelName);
Panel.button(8).Push.Visible = 'on';                                           % Make the "Clear Data" button visible

end


function GetGeoTiff(~,~,panelName)
% Varargin: {'UIControl', 'matlab.ui.eventdata.ActionData', 'panalName'}
global MainWindow ui Geotiff LON LAT

% Geotiff = struct;

% User Interface to select geotiff file
MainWindow.Header.String = 'Select Geotiff file: .tif, .tiff';              % Prompt to select data

layers = MainWindow.Layers;

if isKey(layers, 'DataPoints' )
    vehicle = layers('DataPoints');
    lat = mean(vehicle(:,1));
    lon = mean(vehicle(:,2));
    
    Geotiff = ui.GetGeoTiff(lon, lat );               % Gui to get Geotiff

else
    lat = 44.562692;
    lon = -123.249105;
    Geotiff = ui.GetGeoTiff(lon, lat);                        % Gui to get Geotiff
end

if isempty(Geotiff.Image), return, end                                      % End script if there isn't an input for it to use

% Create geo-referanced mesh for use throughout the program
[m, n, ~] = size(Geotiff.Image);
x = linspace(Geotiff.Data.LongitudeLimits(1), Geotiff.Data.LongitudeLimits(2),n);
y = linspace(Geotiff.Data.LatitudeLimits(2),  Geotiff.Data.LatitudeLimits(1), m);
[LON, LAT] = meshgrid(x,y);

MainWindow.BackGround('LON',LON,'LAT',LAT,'GeoData',Geotiff.Data,'GeoImage',Geotiff.Image);

% Display Geotiff
if any(ismember(keys(MainWindow.Layers),'OcupancyGrid'))
    MakeOcupancyGrid(panelName);
        
end

MainWindow.Header.String = 'Select Data Input or do Color threshold';       % Prompt to select data

% Make Button for Color thresholding visible
Panel = MainWindow.Panels(panelName);
Panel.button(2).Push.Visible = 'on';

end


function GetShapeFile(~,~,panelName)
global MainWindow ui shape

MainWindow.Header.String = 'Select Shape File if Available (.txt) or Cancel Dialog Box';    % Prompt to select data

[fileName, filePath] = ui.getFile('Shape','txt');                           % GUI to get the name and file path of a file
if fileName == 0, return, end                                               % End script if there isn't an input for it to use                                                                % End script if there isn't an input for it to use

shape = ImportWayPoints(filePath,fileName,'\n');                            % Import shape file from .txt file

if shape(1) > 0, shape = fliplr(shape); end                                 % Make sure the waypoints are ordered [Lon Lat], This only works in the north west hemispher

if shape(1,2) ~= shape(end), shape = [shape; shape(1,1), shape(1,2)]; end   % Check if th eshape file is closed
    
MainWindow.Addlayer('Shape', shape, panelName);                             % Dispaly shape on the GUI

Panel = MainWindow.Panels(panelName);                                       % Make buttons usable
Panel.button(5).Push.Visible = 'on';                                        % Clear shape file
Panel.button(5).Push.Enable = 'on';                                         % Enable "clear shape file" button, this will be disabled if the shape file is cleared

end


%% Clear Data
function ClearShapeFile(~,~,panelName)
global MainWindow shape
shape = [];                             % Make shape and empty array
MainWindow.DeleteLayer('Shape');        % Remove layer and its buttons from the GUI
Panel = MainWindow.Panels(panelName);   % Make buttons usable
Panel.button(5).Push.Enable  = 'off';   % Turn off the "Clear shape file" button
end


function ClearAll(~,~)
% Varargin: {'UIControl', 'matlab.ui.eventdata.ActionData', 'panalName'}
disp('Clear all')
run(mfilename)
end


%% Support Functions
function ColorThreshold(~,~)
global MainWindow Geotiff;

colorThresholder(Geotiff.Image);     % Launcth Matlabe's built in color thresholding app
% ColorThresholder;             % Color thresholder that will be included by MATLAB Compiler

MainWindow.Header.String = 'Select Data Input or Make Occupancy Grid';  % Prompt to select data
end


function OpenToolbox(~,event,fncInputts)
% Varargin: {'UIControl', 'matlab.ui.eventdata.ActionData', 'panalName'}
% Execute the selected toolbox function
open = str2func( event.Source.String{ event.Source.Value });
feval(open,fncInputts);

end


function SaveAll(~,~)
% Varargin: {'UIControl', 'matlab.ui.eventdata.ActionData', 'panalName'}
global MainWindow EnvWindow Model DATA

% MainWindow.GUIbuttons( varargin{4},'remove', {'Load Geotiff'});

uisave({'MainWindow','EnvWindow','Model','DATA'},[num2str(yyyymmdd(datetime('today'))),'_EnvironmentalNicheModel']);

disp('Save all -- maybe??')
end




