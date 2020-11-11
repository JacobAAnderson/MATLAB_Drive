% Basic script for reading in StarFish Side Scan Sonar data and displaying it.
% Jacob Anderson
% anderja7@oregonstate.edu
% Dec 12, 2018


clear all
close all
clc

ui = SavedUserInputs(mfilename);                        % Instantiate the Saved User Inputs Class for working with file paths
ui = ui.NewData(false);                                                     % Indicate whether new data should be selescted by the user

% Get Sonar data
% Get file and folder GUI ------------------------------------------------------------------------------------------------
file = ui.GetFile('SonarData','logdoc');                % GUI to get the name and file path of a file

if isempty(file)                                        % Check if the Data Input box was cnaceled
    return                                              % End script if there isn't an input for it to use
end

% Read Logdoc File
Data = Logdoc2mat(file);

% Show Data Header
disp('Data Structure:')
disp(Data.header)

%% Get Data out of structure
% Position data ----------------------------------------------------------------------------------------------------------
utm       = cat(1,Data.position);   % Vehicle location in UTMs [ Easting, Northing, Altitude] 
utmzone   = cat(1,Data.utmZone);    % UTM-zone: Longitude and Latitude bands 
heading   = cat(1,Data.heading);    % Vehicle heading [deg]
speed     = cat(1,Data.speed);      % Vehicle speed   [m/s]
depth     = cat(1,Data.depth);      % Water depth     [m]
portSonar = cat(1,Data.portSonar);  % Port Sonar data structure
starSonar = cat(1,Data.starSonar);  % Starboard Sonar data structure



%% Do stuff with the data
% Display the vehicle's location when the data was collected -------------------------------------------------------------
[lat, lon] = utm2deg(utm(:,1), utm(:,2), utmzone);                          % Conver UTM to lat lon

Geotiff = ui.GetGeoTiff( mean(lon), mean(lat) );                % Gui to get Geotiff

if isempty(Geotiff.Image)                                                        % Check if the Data Input box was cnaceled
    disp('Input box Cancelled')                                             % If no geotiff was selected just plot the location data
    figure('Name','GPS data from Sonar file','NumberTitle', 'off')
    hold on
    scatter(lon,lat,'.')
    p1 = plot(lon(1),lat(1),'*','Color','r');                                % Indicate the starting position
    legend(p1,'Starting Location')
    hold off
    
else
    figure('Name','GPS data from Sonar file','NumberTitle', 'off')          % Plot location data over geotiff
    Geotiff.Show
    hold on 
    scatter(lon,lat,'.')
    p1 = plot(lon(1),lat(1),'*','Color','r');                                % Indicate the starting position
    legend(p1,'Starting Location')
    hold off
end

% Display Sonar image data -----------------------------------------------------------------------------------------------
portEcho = cat(2,portSonar.EchoStrength);                                   % Concatinate the image date from the sonar data structures 
starEcho = cat(2,starSonar.EchoStrength);

image = [flipud(portEcho); starEcho];                                       % Put the port and starboard images together

figure('Name', 'Sonar Image','NumberTitle','off')                           % Show the image
hold on
imshow( image )
l1 = line([0,0],[0,size(image,1)],'color','r');
legend(l1,'Beginning Of Sonar Track', 'Location','northwest')
hold off

Data.image = image;


%% Save file GUI ---------------------------------------------------------------------------------------------------------
% Save Raw Data
disp("\n\nChoose Loaction to Save data as a MATLAB data file '.mat'")

file = ui.saveFile('SideScanSonarData','mat');                              % GUI to get file path to save file
if ~isempty(file)                                                           % Check if the Data Input box was cnaceled
    fprintf('\nSave File Path and name: %s %s \n\n', file);
    save(file,'-struct', 'Data');
end


% Save Image
disp("Choos Location to Save Image as a JPEG")

file = ui.saveFile('image','jpeg');                                         % GUI to get file path to save file
if ~isempty(file)                                                           % Check if the Data Input box was cnaceled
    imwrite(image,file)
end



