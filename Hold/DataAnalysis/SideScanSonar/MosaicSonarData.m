%% Feature Extraction on Side Scan Sonar Images
%  Jacob anderson
%  April 27, 2019

clear all
close all
clc

ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class for working with file paths


%% Get Data sources

newData = true;                                                            % New Data or saved data

if newData
       
    fileType = 'logdoc';                                                    % Type of data file to be processed
    filePath = ui.getDir('SonarData',fileType);                             % GUI to get the filepath to the directory containing the data
    fun      = @(file) Logdoc2mat(file);                                    % Function to be used importing the data
    
    if isempty(filePath)                                                    % Check if the Input box was cnaceled
        return                                                              % End script if there isn't an input for it to use
    end
    
    [geoImage, geoData] = ui.getGeoTiff();                                  % Gui to get Geotiff
    
    Data = BatchDir(filePath,fileType,fun);                                 % Process the directory and return the data
    
    
    file = ui.saveFile('SideScanSonarData','mat');                          % GUI to get file path to save file
    
    if ~isempty(file)                                                       % Check if the Data Input box was cnaceled
        save(file, 'Data');
    end
    
else
    
    file = ui.getFile('SideScanSonarData','mat');                           % GUI to get the name and file path of a file
    
    if isempty(file)                                                        % Check if the Data Input box was cnaceled
        return                                                              % End script if there isn't an input for it to use
    end
    
    [geoImage, geoData] = ui.getGeoTiff();                                  % Gui to get Geotiff
    
    load(file)
end


%% Get Data out of structure and create image
% close all
% clcs

ii = 5;                                                                    % Choose which file to process
% Position data ------------------------------------------------------------------------------------------------------
utm       = Data(ii).position;                                              % Vehicle location in UTMs [ Easting, Northing, Altitude] 
utmzone   = char(Data(ii).utmZone);                                         % UTM-zone: Longitude and Latitude bands 
portSonar = Data(ii).portSonar;                                             % Sonar data structure
starSonar = Data(ii).starSonar;
portEcho  = cat(2,portSonar.EchoStrength);                                  % Concatinate the image date from the sonar data structures 
starEcho  = cat(2,starSonar.EchoStrength);

[lat, lon] = utm2deg(utm(:,1), utm(:,2), utmzone);                          % Conver UTM to lat lon


% Create sonar image --> remove turns
ind = tail(lon,lat);
lon_short = lon(ind);
lat_short = lat(ind);
image = [flipud(portEcho(:,ind)); starEcho(:,ind)];

% Create sonar image --> ignor turns ------
image = [flipud(portEcho); starEcho];   


image(isnan(image)) = 0;                                                    % get rid of NaNs

% Display the vehicle's location when the data was collected ---------------------------------------------------------
figure('Name',sprintf('GPS track: %d',ii),'NumberTitle', 'off')              % Plot GPS coordinates of data
hold on

if ~isempty(geoImage)                                                        % Plot over Geo-tiff is abailalbe
    geoshow(geoImage, geoData)
end

plot(lon,lat)
scatter(lon(1),lat(1),'*')                                                  % Indicate where the path statrted
plot(lon_short,lat_short)
hold off


% Display Sonar image data ----------------------------------------------------------------------------------------------
if mod(ii,2) == 0                                                           % Flip the image around to accomidate the vehicle' direction of travel
    image = fliplr(image);
    image = flipud(image);
end

figure('Name',sprintf('Side Scan Sonar Image --> track: %d',ii),'NumberTitle','off')
hold on
imshow(image)
% plot([1,1],[1,size(image,1)],'r')
hold off

file = ui.saveFile('SideScanSonarImage','mat');                          % GUI to get file path to save file

if ~isempty(file)                                                       % Check if the Data Input box was cnaceled
    filename = Data(ii).header(1);
    filename = filename{1};
    save(file, 'image', 'filename');
    clear filename
end


%% Matlab buitl in Feature Extractors
% close all
% clc

sigma = 6;
image2 = imgaussfilt(image, sigma);

% FAST Features ----------------------------------------------------------------------------
corners = detectFASTFeatures(image2,'MinContrast',0.3);
J = insertMarker(image,corners,'circle');
figure('Name', sprintf('FAST Features --> track: %d, sigma: %d',ii,sigma),'NumberTitle','off')
imshow(J)


% SURF Features ----------------------------------------------------------------------
points_surf = detectSURFFeatures(image2);

figure('Name', sprintf('SURF Features --> track: %d, sigma: %d',ii,sigma),'NumberTitle','off')
imshow( image )
hold on
plot(selectStrongest(points_surf,20))
hold off


% KAZE Features -------------------------------------------------------------------------
points_kaze = detectKAZEFeatures(image2, 'Diffusion','sharpedge',  ...
    'NumScaleLevels', 4, ...
    'Threshold', 0.0000001);  % Find Blobles

figure('Name', sprintf('KAZE Features --> track: %d, sigma: %d',ii,sigma),'NumberTitle','off') 
imshow( image )
hold on
plot(selectStrongest(points_kaze,20))
hold off


% MSER Features -----------------------------------------------------------------------------------
regions = detectMSERFeatures(image2); % Find Regions with different textures

figure('Name', sprintf('MSERF Features --> track: %d, sigma: %d',ii,sigma),'NumberTitle','off')
imshow( image )
hold on
plot(regions,'showPixelList',false,'showEllipses',true);
hold off


% HOG Features ------------------------------------------------------------------------------------
[featureVector,hogVisualization] = extractHOGFeatures(image2);

figure('Name', sprintf('HOG Features --> track: %d, sigma: %d',ii,sigma),'NumberTitle','off')
imshow( image )
hold on
plot(hogVisualization);
hold off


% Harris Features ----------------------------------------------------------------------------------
points_harris = detectHarrisFeatures(image2);

figure('Name', sprintf('Harris Features --> track: %d, sigma: %d',ii,sigma),'NumberTitle','off') 
imshow( image )
hold on
plot(selectStrongest(points_harris,20))
hold off


%% Plot data track
utm       = cat(1,Data.position);                                              % Vehicle location in UTMs [ Easting, Northing, Altitude] 
utmzone   = char(cat(1,Data.utmZone));                                         % UTM-zone: Longitude and Latitude bands 

[lat, lon] = utm2deg(utm(:,1), utm(:,2), utmzone);                          % Conver UTM to lat lon

figure('Name','GPS track','NumberTitle', 'off')              % Plot GPS coordinates of data
hold on

if ~isempty(geoImage)                                                        % Plot over Geo-tiff is abailalbe
    geoshow(geoImage, geoData)
end

plot(lon,lat)
scatter(lon(1),lat(1),'*')                                                  % Indicate where the path statrted
hold off


%% Custom Featur Extraction 
domain = [-pi/1, pi/1];
sigma = 5;
kernalSize = [10, 10];                                                     % [row, column] --> [y , x ]

% Make the Kernal ----------------------------------------------------------------------------------------------------
kernal = @(x,y) exp( -(x.*x + y.*y)./(sigma * sigma) );                   % Kernal Function
%kernal = @(x,y) ( cos( sqrt(x.*x + y.*y) ) );                                                  % Kernal Function

q = integral2(kernal,domain(1),domain(2), domain(1),domain(2));             % Normalize the kernal function
kernal = @(x,y) ( kernal(x,y) - q/4 );                                      % Normalized Kernal function

xx = linspace( domain(1), domain(2), kernalSize(2));                        % Window Size
yy = linspace( domain(1), domain(2), kernalSize(1));

[xx, yy] = meshgrid(xx,yy); 

KK = kernal(xx,yy);                                                         % Kernal as a matrix for convolution
KK = KK ./ numel(KK);

% Do Convolution and post processing ----------------------------------------------------------------------------------
% C = conv2(image,KK);                                                        % Convolution
C = imgaussfilt(image, 15);
% C = C./ max(C);

figure('Name','convolved Imange')
imshow(C)

imSize = size(image);                                                       % Trim the convolved image back down to the origonal size
cSize  = size(C);
trim = floor( (cSize - imSize)./2 );

C(imSize(1) - trim(1): imSize(1), : ) = [];
C(:, imSize(2) - trim(2): imSize(2) ) = [];
C(1:trim(1), : ) = [];
C(: , 1:trim(2)) = [];


C(C<15000) = NaN;

% Visulaize results ----------------------------------------------------------------------------------------------------
figure('Name','Sonar Image','NumberTitle','off')
imshow(image)
hold on
x = [0, 0, kernalSize(2), kernalSize(2), 0];
y = [0, kernalSize(1), kernalSize(1), 0, 0];
plot(x,y, 'r-', 'LineWidth', 1);
hold off

figure('Name','Kernal','NumberTitle','off')
mesh(KK)
xlabel('X')
ylabel('Y')

figure('Name','Image Convolution Responce','NumberTitle','off')
mesh(C)
view(30,45)

 
 
%%  Find Sonar traughff
close all
clc

KK = [ 2,  4,  4,  4,  4,  2;
       2,  4,  4,  4,  4,  2;
       4,  6,  8,  8,  6,  4;
       0,  0,  0,  0,  0,  0;
      -4, -6, -8, -8, -6, -4;
      -2, -4, -4, -4, -4, -2;
      -2, -4, -4, -4, -4, -2];

KK = [KK,KK];
  
%KK =flipud(KK);  
  
C = conv2(image,KK);

A = zeros(size(C));


A(C >  max(C(:)) * 0.4 ) = 1;
A(C <  min(C(:)) * 0.2 ) = 1;

[row,~] = size(A);

A(1:10,:) = 0;
A(row/2 -20 : row/2 + 20, :) = 0;
A(row - 10: row,:) = 0;

imshow(A)

pause(2)

% KK = [ 0,  0,  0,  0,  0,  0;
%       -1, -1, -1, -1, -1, -1;
%        1,  1,  1,  1,  1,  1;
%       -1, -1, -1, -1, -1, -1;
%       -0,  0,  0,  0,  0,  0];
  
C = conv2(A,KK);

imshow(C)

% h = fspecial('average',[5, 5]);                            % Create an averaging filter to smoth edges
% A = imfilter(A,h);                                     % Applie the filter

% K = [ 2,  2,  2,  2;
%       4,  4,  4,  4;
%       6,  6,  6,  6;
%       0,  0,  0,  0;
%      -6, -6, -6, -6;
%      -4, -4, -4, -4;
%      -2, -2, -2, -2];
%  
% KK = [K,K,K,K,K,K,K,K,K];   
% A = conv2(A,KK);

% imshow(A)

% B = bwboundaries(A);
% 
% Len = cellfun('length', B);
% [Len, I] = sort(Len, 'descend' );
% 
% imshow(image)
% hold on
% 
% boundary = B{I(1)};
% plot(boundary(:,2), boundary(:,1), 'LineWidth', 2)
% 
% boundary = B{I(2)};
% plot(boundary(:,2), boundary(:,1), 'LineWidth', 2)


%% Functions
function ind = tail(lon,lat)
c = polyfit(lon,lat,1);
lat_est = polyval(c,lon);

pathDev = lat-lat_est;
mu = mean(pathDev);
si = std(pathDev);

out = find(pathDev > mu+si | pathDev < mu-si);

if any(out(:) == 1 )
    ind = max(out):numel(lon);
else
    ind = 1: min(out);
end

end
