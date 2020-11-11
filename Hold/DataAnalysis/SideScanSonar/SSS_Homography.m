

close all
clear all
clc

ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class

% Get file with file path GUI ----------------------------------------------------------------------------------------------
file1 = ui.getFile('SSS_Data1','mat');                                      % GUI to get the name and file path of a file

if isempty(file1)                                                            % Check if the Data Input box was cnaceled
    return                                                                  % End script if there isn't an input for it to use
end

file2 = ui.getFile('SSS_Data2','mat');                                      % GUI to get the name and file path of a file

if isempty(file1)                                                            % Check if the Data Input box was cnaceled
    return                                                                  % End script if there isn't an input for it to use
end


% Load Data

Data1 = load(file1);

Data2 = load(file2);

clear file1 file2

%% Get Data out of structure
% Data file 1
utm1       = cat(1,Data1.position);   % Vehicle location in UTMs [ Easting, Northing, Altitude] 
utmzone1   = cat(1,Data1.utmZone);    % UTM-zone: Longitude and Latitude bands 
heading1   = cat(1,Data1.heading);    % Vehicle heading [deg]
% speed1     = cat(1,Data1.speed);      % Vehicle speed   [m/s]
depth1     = cat(1,Data1.depth);      % Water depth     [m]
portSonar1 = cat(1,Data1.portSonar);  % Port Sonar data structure
starSonar1 = cat(1,Data1.starSonar);  % Starboard Sonar data structure

% Data file 2
utm2       = cat(1,Data2.position);   % Vehicle location in UTMs [ Easting, Northing, Altitude] 
utmzone2   = cat(1,Data2.utmZone);    % UTM-zone: Longitude and Latitude bands 
heading2   = cat(1,Data2.heading);    % Vehicle heading [deg]
% speed2     = cat(1,Data2.speed);      % Vehicle speed   [m/s]
depth2     = cat(1,Data2.depth);      % Water depth     [m]
portSonar2 = cat(1,Data2.portSonar);  % Port Sonar data structure
starSonar2 = cat(1,Data2.starSonar);  % Starboard Sonar data structure



%% Do stuff with the data
% Display the vehicle's location when the data was collected -------------------------------------------------------------
[lat1, lon1] = utm2deg(utm1(:,1), utm1(:,2), utmzone1);                          % Conver UTM to lat lon
[lat2, lon2] = utm2deg(utm2(:,1), utm2(:,2), utmzone2);                          % Conver UTM to lat lon

[geoImage, geoData] = ui.getGeoTiff( mean(lon1), mean(lat1) );                % Gui to get Geotiff

if isempty(geoImage)                                                        % Check if the Data Input box was cnaceled
    disp('Input box Cancelled')                                             % If no geotiff was selected just plot the location data
    figure('Name','GPS data from Sonar file','NumberTitle', 'off')
    hold on
    scatter(lon1,lat1,'.','g')
    scatter(lon2,lat2,'.','c')
    p1 = plot(lon1(1),lat1(1),'*','Color','r');                                % Indicate the starting position
    plot(lon2(1),lat2(1),'*','Color','r');
    legend(p1,'Starting Location')
    hold off
    
else
    figure('Name','GPS data from Sonar file','NumberTitle', 'off')          % Plot location data over geotiff
    geoshow(geoImage, geoData)
    hold on 
    scatter(lon1,lat1,'.','g')
    scatter(lon2,lat2,'.','c')
    plot(lon2(1),lat2(1),'*','Color','r');
    p1 = plot(lon1(1),lat1(1),'*','Color','r');                                % Indicate the starting position
    legend(p1,'Starting Location')
    hold off
end

% Display Sonar image data -----------------------------------------------------------------------------------------------
portEcho1 = cat(2,portSonar1.EchoStrength);                                   % Concatinate the image date from the sonar data structures 
starEcho1 = cat(2,starSonar1.EchoStrength);

image1 = [flipud(portEcho1); starEcho1];


portEcho2 = cat(2,portSonar2.EchoStrength);                                   % Concatinate the image date from the sonar data structures 
starEcho2 = cat(2,starSonar2.EchoStrength);

image2 = [flipud(portEcho2); starEcho2];              

figure('Name','SSS Images','NumberTitle','off')
subplot(2,1,1)
imshow(image1)
subplot(2,1,2)
imshow(image2)



%% Create Homography Matrix
t1 = [mean(utm1(:,1)), mean(utm1(:,2)), 0];     % Image 1 Translation 
t2 = [mean(utm2(:,1)), mean(utm2(:,2)), 0];     % Image 2 Translation

yaw1 = mean(heading1);                          % Image 1 orientation
yaw2 = mean(heading2);                          % Image 2 orientation

d = mean(depth2);
n = [0,0,1];

Rz = @(yaw) [cosd(yaw) -sind(yaw) 0; sind(yaw) cosd(yaw) 0; 0 0 1];

R1 = Rz(yaw1);
R2 = Rz(yaw2);

I = eye(3);

H = R1*( I - (t2-t1)*n'/d) * R2';

disp("")
disp('Homography Matrix:')
disp(H)

% Save Homograpy matrix in xml format
xml1 = sprintf('<?xml version="1.0"?>\n<opencv_storage>\n<H13 type_id="opencv-matrix">\n  <rows>3</rows>\n  <cols>3</cols>\n  <dt>d</dt>\n  <data>');
hom  = sprintf('\n\t%f  %f  %f',H);
xml2 = sprintf('</data></H13>\n</opencv_storage>');

file = ui.saveFile('homography_matrix','xml');                                           % GUI to get file path to save file

if isempty(file)                                                            % Check if the Data Input box was cnaceled
    return                                                                  % End script if there isn't an input for it to use
end

fileID = fopen(file,'w');
fprintf(fileID, '%s%s%s', xml1, hom, xml2);
fclose(fileID);

