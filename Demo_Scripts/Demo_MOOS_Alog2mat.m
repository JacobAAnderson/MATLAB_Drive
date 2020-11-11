clear all
close all
clc

UserInputs = SavedUserInputs(mfilename);                                    % Instantiate the Saved User Inputs Class

% Get file and folder GUI ----------------------------------------------------------------------------------------------
[fileName, filePath] = UserInputs.getFile('alog','alog');                   % GUI to get the name and file path of a file
fprintf('\nReading In:\n%s %s \n\n', filePath, fileName);

if fileName == 0                                                            % Check if the Data Input box was cnaceled
    disp('Input box Cancelled')
    return                                                                  % End script if there isn't an input for it to use
end

Data = MOOS_Alog2mat(fileName, filePath);                                   % Import MOOS Alog

clear fileName filePath


%% Do usefull things with data
% Find get lat - lon values -------------------------------------------------------------------------------------------------------------
close all

gps = cat(1,Data.gps);
%nav = cat(1,Data.nav);
alt = cat(1,Data.alt);


% Plot data times
% figure('Name','Times That data was Collected','NumberTitle','off')
% hold on
%
% times = nav(:,1);
% times = cat(1,times{:});
% plot(times, zeros(size(times)), '.','Color', 'b')
%
% times = gps(:,1);
% times = cat(1,times{:});
% plot(times, ones(size(times))* 0.01, '*','Color','r')
%
% times = alt(:,1);
% times = cat(1,times{:});
% plot(times, ones(size(times))* 0.02, '*','Color','g')
%
% ylim([-0.5 0.5])
%
% hold off



%% Get Geotiff of the area
lon = gps(:,2);                                                             % Find center of the map
lon = cat(1,lon{:});
lon0 = mean(lon(:));

lat = gps(:,3);
lat = cat(1,lat{:});
lat0 = mean(lat(:));

[geoImage, geoData] = UserInputs.getGeoTiff(lon0,lat0);


% Plot Geotiff and GPS data ----------------------------------------------------
figure('Name','Geo-tiff','NumberTitle','off')
hold on
geoshow(geoImage, geoData)
plot(lon0, lat0, '*', 'Color', 'r')
plot(lon,  lat,  '*', 'Color', 'g')
hold off






