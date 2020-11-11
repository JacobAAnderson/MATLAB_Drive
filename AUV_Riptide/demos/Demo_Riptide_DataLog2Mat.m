% Demo RiptideDataLog2Mat


close all
clear all
clc


ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(true);                                                      % Indicate whether new data should be selescted by the user

% Get file with file path GUI ----------------------------------------------------------------------------------------------
data_file = ui.GetFile('Data','csv');                                       % Get file path GUI
if isempty(data_file), return, end                                          % End script if there isn't an input for it to use
 

offset_file = ui.GetFile('Off Set','txt');                                  % Get file path GUI

% ui = ui.NewData(false); 
geotiff = ui.GetGeoTiff;
if isempty(geotiff.Image), return, end                                      % End script if there isn't an input for it to use


%% Extract data
filters = {'GPS_fix', 0;                                                    % Data filter {"Data_Field", value to be filterred out}
           'ALT_ALTITUDE', 0};  

Data = Riptide_DataLog2Mat(data_file, offset_file, filters); %, filters);

rt_lat = Data.vehicle(:,1);
rt_lon = Data.vehicle(:,2);

bathy_lat = Data.bathymetry(:,1);
bathy_lon = Data.bathymetry(:,2);

mission = Data.mission;

idel = strcmp(mission, 'idel');


%% Show some data

geotiff.Show
hold on
plot(rt_lon,rt_lat,'.b')
%plot(bathy_lon,bathy_lat,'.w')
hold off
