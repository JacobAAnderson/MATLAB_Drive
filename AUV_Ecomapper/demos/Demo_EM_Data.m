% Demo EM_Data Class

close all
clear all
clc

ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(true);                                                      % Indicate whether new data should be selescted by the user

file = ui.GetFile('Ecomapper','log');                                       % GUI to get the name and file path of a file
if isempty(file), return, end                                               % End script if there isn't an input for it to use

Geotiff = ui.GetGeoTiff;                                                    % Get Geotiff GUI

emData = EM_Data(file);                                                     % Instanciate the class

emData = emData.FilterAltimiter(2);                                         % Low Pass Filter for Raw Altimiter Data, Cut data points that excced 2 std
emData = emData.FilterHeading(1);                                           % Low Pass Filter for Raw Compass Data, Cut data points that excced 1 std

emData = emData.GetVehiclePaths;                                            % Generate the vehicle's path in lat-lon, utm and dead reckoning
emData = emData.Altimeter_Path;                                             % Generate the Bathymetry readings in lat-lon, utm, dead reckoning

station_ID = 9410079;                                                       % NOAA station ID for Catlina IS
emData = emData.TideCorrection(station_ID);                                 % Get the Tide Level that correspondes to the when the data was colected 
tide_Height = emData.RawData.tide;                                          % Access Tide data. This has not been applied to the Raw Data

emData.PlotPaths(Geotiff)                                                   % Plot Pathes over the geotiff

clear file station_ID


%%
clc
close all

lat_lon = emData.RawData.vehicle(:,1:2);
raw_alt = emData.RawData.attitude(:,4);

for k = 1: 30

[idx,C] = kmeans(lat_lon,k);

alt(k) = 0;

figure

for ii = 1:k
   
    alt(ii) = std(raw_alt(idx == ii));

    plot(lat_lon(idx == ii, 1), lat_lon(idx == ii, 2), '.')
    hold on

end


fprintf('\nK: %d\tAltimiter STD: %f\n',k, mean(alt))

end




