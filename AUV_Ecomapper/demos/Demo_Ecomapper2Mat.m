% Demo Script for the Ecomapper2Mat Function
close all
clear all
clc

% Get get data sources ----------------------------------------------------
ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class

file = ui.getFile('EcomapperData','log');                                   % GUI to get the name and file path of a file

if isempty(file)                                                            % Check if the Data Input box was cnaceled
    return                                                                  % End script if there isn't an input for it to use
end

% Import Data -------------------------------------------------------------
Data = EcoMapperLog2Mat(file,'dvl filter');

disp('Data Structure:')
disp(Data.header)


%% Get Data out of data structure
vehicle         = cat(1, Data.vehicle);
bathymetry      = cat(1, Data.bathymetry);
waterParameters = cat(1, Data.wqData);

lon   = vehicle(:,2);
lat   = vehicle(:,1);
depth = vehicle(:,3);
temp  = waterParameters(:,1);


%% Do Stuff with Data

% Show GPS track ----------------------------------------------------------
[geoImage, geoData] = ui.getGeoTiff(lon, lat);                                        % Gui to get Geotiff ---> (lon, lat) is optional for locating file on your computer but needed for retriving map from web map server

figure('Name','GPS track','NumberTitle','off')
hold on
if ~isempty(geoImage)                                                        % Check if the Data Input box was cnaceled
    geoshow(geoImage, geoData)
end
plot(lon,lat,'g')
hold off
% export_fig missionPath.png    

% Plot a water parameter --------------------------------------------------
figure('Name','Water Temperature','NumberTitle','off')
scatter3( lon, lat, -depth, '.', 'CData', temp )
colorbar