
% 3D Reconstruction from Sonar Data


% clear all
close all
clc

UserInputs = SavedUserInputs(mfilename);                % Instantiate the Saved User Inputs Class for working with file paths

% Get Sonar data ----------------------------------------------------------------------------------------------
file = UserInputs.getFile('SonarData','logdoc');        % GUI to get the name and file path of a file
fprintf('\nFile Path and name: %s %s \n\n', file);

if isempty(file)                                        % Check if the Data Input box was cnaceled
    return                                              % End script if there isn't an input for it to use
end

Data = Logdoc2mat(file);                                % Get Data


% % Save data?? --> Uncoment to save the data -------------------------------------------------------------------------------------
% file = UserInputs.saveFile('SideScanSonarData','mat');          % GUI to get file path to save file
%
% if ~isempty(file)                                               % Check if the Data Input box was cnaceled
%     fprintf('\nSave File Path and name: %s %s \n\n', file);
%     save(file, 'Data');
% end


% Get Data out of structure ------------------------------------------------------------------------------------------------------
utm       = cat(1,Data.position);           % Vehicle location in UTMs [ Easting, Northing, Altitude]
utmzone   = cat(1,Data.utmZone);            % UTM-zone: Longitude and Latitude bands
heading   = cat(1,Data.heading);            % Vehicle heading in degrees
speed     = cat(1,Data.speed);              % Vehicle speed in meters per second
depth     = cat(1,Data.depth);              % Water depth in meters
portSonar = cat(1,Data.portSonar);          % Sonar data structure
starSonar = cat(1,Data.starSonar);

portEcho  = cat(2,portSonar.EchoStrength);  % Concatinate the image date from the sonar data structures
starEcho  = cat(2,starSonar.EchoStrength);


image = [flipud(portEcho); starEcho];       % Put the port and starboard images together


%% Display Sonar image data -----------------------------------------------------------------------------------------------
figure('Name', 'Sonar Image','NumberTitle','off')
hold on
imshow( image )
l = line([1,1],[1,size(image,1)],'Color','r');      % Line to indicate the which side the vehicle started on
legend(l,'Sart', 'Location','southwest')
hold off
clear l


%% Display the vehicle's location when the data was collected -------------------------------------------------------------
[lat, lon] = utm2deg(utm(:,1), utm(:,2), utmzone);                          % Conver UTM to lat lon
% 
% [geoImage, geoData] = UserInputs.getGeoTiff( mean(lon), mean(lat) );        % Gui to get Geotiff
% 
% if isempty(geoImage)                                                        % Check if the Data Input box was cnaceled
%     disp('Input box Cancelled')                                             % If no geotiff was selected just plot the location data
%     figure('Name','GPS data from Sonar file','NumberTitle', 'off')
%     hold on
%     scatter(lon,lat,'.')
%     plot(lon(1),lat(1),'*','Color','r')                                     % Indicate the starting position
%     hold off
% else
%     figure('Name','GPS data from Sonar file','NumberTitle', 'off')          % Plot location data over geotiff
%     geoshow(geoImage, geoData)
%     hold on
%     scatter(lon,lat,'.')
%     plot(lon(1),lat(1),'*','Color','r')
%     hold off
% end


%% Get trench
portTrench = FindTrench(portEcho);
starTrench = FindTrench(starEcho);

figure('Name','Sonar Image with Trench','Numbertitle','off')
subplot(2,1,1)
imshow(portEcho)
hold on
plot(portTrench, '.', 'Color', 'r')
hold off
drawnow

subplot(2,1,2)
imshow(starEcho)
hold on
plot(starTrench, '.', 'Color', 'r')
hold off
drawnow



%% Shape from shade
maxRange = 30;

deep = -10;

trench(1,:) = portTrench;
trench(2,:) = starTrench;

s = vdist(lat(1),lon(1),lat(end),lon(end));

[zPort,zStar] = BuildHeightMap(image, deep, utm, s, maxRange, xOffset, trench);






%% Find Trench
function trench = FindTrench(echoImage)

echoImage(isnan(echoImage)) = 0;

[row, column] = size(echoImage);
trench = zeros(1,column);


filter    = fspecial('average',[20,20]);
echoImage = imfilter(echoImage,filter);

filter    = fspecial('disk',20);
echoImage = imfilter(echoImage,filter);

mu    = nanmean(echoImage(:));
% sigma = nanstd(echoImage(:));


% Find the edge of the sonar trough ------------------------------------------------------------------------------------
for jj = 3: column - 2
    
    ii = 10;
    while true
        m = nanmean(echoImage( ii : ii+10, jj-2:jj+2 ));
        
        if m > mu %+ sigma 
            
            if jj == 3
                trench(1) = ii;
                trench(2) = ii;
                trench(3) = ii;

            else
                trench(jj) = 0.7*trench(jj-1) + 0.3*ii;
            end
            
            break
            
            
        elseif ii >= row - 101
            trench(jj) = trench(jj-1);
            break
            
        else
            ii = ii+1;
            
        end
    end
    
end

trench(column)   = ii;
trench(column-1) = ii;

% Smooth out the demarcation line -------------------------------------------------------------------------------------------------------------
for jj = 1:5
    for ii = (column-1): -1: 1
        trench(ii) = trench(ii)*0.2 + 0.8*trench(ii+1);
    end
    
    for ii = 2: column
        trench(ii) = 0.7*trench(ii-1) + 0.3*trench(ii);
    end
end

trench = ceil(trench);

end


