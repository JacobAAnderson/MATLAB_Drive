
clc
clear all
close all

UserInputs = SavedUserInputs(mfilename);                                    % Instantiate the Saved User Inputs Class

%% Get Sonar data
% Get file and folder GUI ----------------------------------------------------------------------------------------------
[fileName, filePath] = UserInputs.getFile('SonarData','logdoc');            % GUI to get the name and file path of a file
fprintf('\nFile Path and name: %s %s \n\n', filePath, fileName);

if fileName == 0                                                            % Check if the Data Input box was cnaceled
    return                                                                  % End script if there isn't an input for it to use
end

file = fullfile(filePath, fileName);

Data = Logdoc2mat(file);


% Save file GUI ---------------------------------------------------------------------------------------------------------
% [fileName, filePath] = UserInputs.saveFile('SideScanSonarData','mat');                   % GUI to get file path to save file
%
% if fileName ~= 0                                                            % Check if the Data Input box was cnaceled
%     fprintf('\nSave File Path and name: %s %s \n\n', filePath, fileName);
%     file = fullfile(filePath, fileName);
%     save(file, 'Data');
% end

clear fileName filePath file


%% Get Data out of structure
% Position data ------------------------------------------------------------------------------------------------------
position = cat(1,Data.position);
XX = position(:,1);
YY = position(:,2);

utmzone  = char(cat(1,Data.utmZone));

[lat, lon] = utm2deg(XX, YY, utmzone);    % Conver UTM to lat lon

% plot on a Geotiff --------------------
[geoImage, geoData] = UserInputs.getGeoTiff( mean(lon), mean(lat) );        % Gui to get Geotiff


% Heading data --------------------------------------------------------------------------------------------------------
heading  = cat(1,Data.heading);
depth    = cat(1,Data.depth);

% speed = cat(1,Data.speed);

% Sonar data
portSonar    = cat(1,Data.portSonar);
starSonar    = cat(1,Data.starSonar);


%% Mosaic Sonar data
plotGeo = true();

% Creat mamp for the mosaic
portResolution = double( cat(1,portSonar.Resolution) );
starResolution = double( cat(1,starSonar.Resolution) );

sonarResolution = min( min(portResolution(:),starResolution(:)) );
sonarSwath = sonarResolution * length(portSonar(1).EchoStrength);

padding = 1.7 * sonarSwath;

origon   = [min( position(:,1)) - padding, min( position(:,2)) - padding];
topright = [max( position(:,1)) + padding, max( position(:,2)) + padding];

gridsize = floor( fliplr( 2 * (topright - origon)./sonarResolution ));

map = NaN(gridsize);

clear portResolution starResolution gridsize

utmzone = utmzone(1,:);

% Process Data
figure('Name','Sonar Coverage','NumberTitle','off')

if ~isempty(geoImage) && plotGeo
    geoshow(geoImage, geoData)
end

hold on

for ii = 1: length(Data)
    
    % Get the locations sonar return in world coordinate frame
    [xy_Port, xy_Star] = SonarReturn_XY( XX(ii), YY(ii), heading(ii), depth(ii), sonarResolution, length(portSonar(ii).EchoStrength) );
    
    index = GetIndecies(xy_Port, map, origon, topright); % Get the map indecies cooresponding to the port side sonar readings
    map(index) = portSonar(ii).EchoStrength;
    
    index = GetIndecies(xy_Star, map, origon, topright); % Get the map indecies cooresponding to the starboard side sonar readings
    map(index) = starSonar(ii).EchoStrength;
    
    
    % Plot location of data points --------------------------------------------------------------------------------------------------
    if plotGeo
        [Vlat, Vlon] = utm2deg(XX(ii), YY(ii), utmzone);
        plot(Vlon, Vlat, '.', 'Color','b')
        
        zone = repmat(utmzone,size(xy_Star(:,1)));
        [Plat, Plon] = utm2deg(xy_Port(:,1),xy_Port(:,2),zone);
        plot(Plon, Plat,'.', 'color', [ 0.9294    0.6902    0.1294 ])
        
        [Slat, Slon] = utm2deg(xy_Star(:,1),xy_Star(:,2),zone);
        plot(Slon, Slat, '.', 'Color', [ 0.4706    0.6706    0.1882 ])
    end
    
end

hold off

figure
imshow(flipud(map),'DisplayRange',[0 1])


%% Ride'n dirty

% close all
portEcho = cat(2,portSonar.EchoStrength);
starEcho = cat(2,starSonar.EchoStrength);

map2 = [flipud(portEcho); starEcho];

figure
imshow( map2 )



%% Functions ================================================================================================================================
function index = GetIndecies(xy, matrix, origon, topright)

[r, c] = size(matrix);

res = (topright - origon) ./ [c, r];                % Size of array elemant

column = floor( (xy(:,1) - origon(1)) / res(1) );   % X position
row    = floor( (xy(:,2) - origon(2)) / res(2) );   % Y position

index = sub2ind([r, c], row, column );              % Get liner index for multiple assignmants later

end


function [xy_Port, xy_Star] = SonarReturn_XY(x, y, theta, depth, sonarRes, numbReadings)

L = depth*cosd(30);         % horizontal distance from underneath the vehicle to the frist sonar reading
l = (1:1: numbReadings).* double(sonarRes) + L;      % Vector representing the sonar readings

xy_Port = [ x + l'*sind(theta - 90.0), y + l'*cosd(theta - 90.0)];

xy_Star = [x + l'*sind(theta + 90.0), y + l'*cosd(theta + 90.0)];

end


