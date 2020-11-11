% Demo Script for the Ecomapper2Mat Function
close all
clear all
clc

close all
clear all
clc

ui = SavedUserInputs(mfilename);                                    % Instantiate the Saved User Inputs Class

fileType = 'log';
fun = @(file) EcoMapperLog2Mat(file,'dvl filter');

% Get folder GUI ------------------------------------------------------------------------------------------------------
filePath = ui.getDir('Data',fileType);                                 % Get folder GUI

if isempty(filePath)                                                        % Check if the Input box was cnaceled
    disp('Get Directory Input box Cancelled')
    return                                                                  % End script if there isn't an input for it to use
else
    fprintf('\nFile Path to the Dir: %s \n\n', filePath);    
end

Data = BatchDir(filePath,fileType,fun);

disp('Data Structure:')
disp(Data(1).header)


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
% Enlarge the viewing area
dLon = max(lon) - min(lon);
dLat = max(lat) - min(lat);

tiff_lon = [min(lon) - 0.2 * dLon, max(lon) + 0.2 * dLon];
tiff_lat = [min(lat) - 0.2 * dLat, max(lat) + 0.2 * dLat];

% Get Tiff
[geoImage, geoData] = ui.getGeoTiff(tiff_lon, tiff_lat);                                        % Gui to get Geotiff ---> (lon, lat) is optional for locating file on your computer but needed for retriving map from web map server

figure('Name','GPS track','NumberTitle','off')
hold on
if ~isempty(geoImage)                                                        % Check if the Data Input box was cnaceled
    geoshow(geoImage, geoData)
end
plot(lon,lat,'g')
hold off
    

% Plot a water parameter --------------------------------------------------
figure('Name','Water Parameter','NumberTitle','off')
scatter3( lon, lat, -depth, '.', 'CData', temp )
colorbar