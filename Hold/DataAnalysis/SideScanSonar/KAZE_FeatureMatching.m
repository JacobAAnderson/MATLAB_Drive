% Feature Matching
% Jacob Anderson
% May 7, 2019

close all
clear all
clc

ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class


%% Get Data
file1 = ui.getFile('Data1','mat');                                          % GUI to get the name and file path of a file

if isempty(file1)                                                           % Check if the Data Input box was cnaceled
    return                                                                  % End script if there isn't an input for it to use
end


file2 = ui.getFile('Data2','mat');                                          % GUI to get the name and file path of a file

if isempty(file2)                                                           % Check if the Data Input box was cnaceled
    return                                                                  % End script if there isn't an input for it to use
end


[geoImage, geoData] = ui.getGeoTiff;                                        % Gui to get Geotiff ---> (lon, lat) is optional for locating file on your computer but needed for retriving map from Google / web map server

if isempty(geoImage)                                                        % Check if the Data Input box was cnaceled
    return                                                                  % End script if there isn't an input for it to use
end


data1    = load(file1);                                                     % Load data from file 1
image1   = data1.image;                                                     % Side Scan Soanr Image
utm1     = cat(1,data1.position);                                           % Vehicle location in UTMs [ Easting, Northing, Altitude] 
utmzone1 = cat(1,data1.utmZone);                                            % UTM-zone: Longitude and Latitude bands 
heading1 = cat(1,data1.heading);                                            % Vehicle heading [deg]
depth1   = cat(1,data1.depth);                                              % Water depth     [m]

data2    = load(file2);                                                     % Load data from file 2
image2   = data2.image;                                                     % Side Scan Soanr Image
utm2     = cat(1,data2.position);                                           % Vehicle location in UTMs [ Easting, Northing, Altitude] 
utmzone2 = cat(1,data2.utmZone);                                            % UTM-zone: Longitude and Latitude bands 
heading2 = cat(1,data2.heading);                                            % Vehicle heading [deg]
depth2   = cat(1,data2.depth);                                              % Water depth     [m]

clear file1 file2


%% Find Features
points1 = detectKAZEFeatures(image1, 'Diffusion','region', 'NumOctaves', 4); %, 'NumScaleLevels', 4, 'Threshold', 0.0000001);  % Find Blobles
points2 = detectKAZEFeatures(image2, 'Diffusion','region', 'NumOctaves', 4); %, 'NumScaleLevels', 4, 'Threshold', 0.0000001); 


%% Reduce Features

% Eliminate weak features -------------------------------------------------
points1_filtered = selectStrongest(points1, 30);
points2_filtered = selectStrongest(points2, 30);

% Eliminate features in sonar traugh --------------------------------------
midLine = floor(size(image1,1)/2);
points1_filtered = ClearTrench(points1_filtered, midLine);
points2_filtered = ClearTrench(points2_filtered, midLine);

% points1_filtered = MergPoints(points1_filtered);
% points2_filtered = MergPoints(points2_filtered);


%% Dispaly features
figure('Name', "Features 1",'NumberTitle','off')
%subplot(2,1,1)
imshow( image1 )
hold on
plot(points1_filtered)
hold off
export_fig imag1.png

figure('Name', "Features 2",'NumberTitle','off')
%subplot(2,1,2)
imshow( image2 )
hold on
plot(points2_filtered)
hold off
export_fig image2.png


%% Get feature Coordinates

xy1 = GetFeatureCoords( cat(1,points1_filtered.Location), depth1, heading1, midLine, utm1);
xy2 = GetFeatureCoords( cat(1,points2_filtered.Location), depth2, heading2, midLine, utm2);

% xy1 = round( cat(1,points1_filtered.Location) );                            % X,Y coordinates of feature on SSS image [pixels]
% xy2 = round( cat(1,points2_filtered.Location) );
% 
% height1 = depth1( xy1(:,1) );                                               % Height of the vehicle above sea floor [meters]
% height2 = depth1( xy2(:,1) );
% 
% theta1 = heading1( xy1(:,1) );                                              % Heading of the vehicle coorisponding to the features [deg]
% theta2 = heading2( xy2(:,1) );
% 
% dist1 = (xy1(:,2) - midLine) .* 30/midLine;                                 % Range [meters]
% dist2 = (xy2(:,2) - midLine) .* 30/midLine;
% 
% dist1 = sqrt(dist1.*dist1 - height1.*height1) .* sign(dist1);               % Horizontal distance of the feature [meters]
% dist2 = sqrt(dist2.*dist2 - height2.*height2) .* sign(dist2);
% 
% xy1 = utm1( xy1(:,1), : );                                                  % Location of the vehicle [utms "meters"]
% xy2 = utm1( xy2(:,1), : );
% 
% theta1 = theta1 + 90;                                                       % Orientation of the feature [deg]
% theta2 = theta2 + 90;
% 
% theta1( theta1 < 0 )   = theta1( theta1 < 0 ) + 360;                        % Keep  360 >= theta >= 0
% theta1( theta1 > 360 ) = theta1( theta1 > 360 ) - 360;
% 
% theta2( theta2 < 0 )   = theta2( theta2 < 0 ) + 360;
% theta2( theta2 > 360 ) = theta2( theta2 > 360 ) - 360;
% 
% xy1(:,1) = xy1(:,1) + dist1 .* sind( theta1 );                              % Global position of the features
% xy1(:,2) = xy1(:,2) + dist1 .* cosd( theta1 );
% 
% xy2(:,1) = xy2(:,1) + dist2 .* sind( theta2 ); 
% xy2(:,2) = xy2(:,2) + dist2 .* cosd( theta2 );


%% Plot the vehicles path
utmzone = utmzone1(1,:);

figure('Name','Sonar Track1 with Landmarks','Numbertitle','off')
geoshow(geoImage, geoData)
hold on

[lat, lon] = utm2deg( utm1(:,1), utm1(:,2), utmzone1 );                     % Conver UTM to lat lon
p1 = plot(lon,lat,'y');                                                     % Plot Vehicle Path
scatter(lon(1), lat(1),'d','Cdata', [0,0,0])

[lat, lon] = utm2deg( xy1(:,1), xy1(:,2), repmat(utmzone,size(xy1,1),1) );  % Conver UTM to lat lon
s1 = scatter(lon, lat,'+','cdata',[1  0  1]);                                     % Plot Landmarks

hold off


figure('Name','Sonar Track2 with Landmarks','Numbertitle','off')
geoshow(geoImage, geoData)
hold on

[lat, lon] = utm2deg(utm2(:,1), utm2(:,2), utmzone2);                       % Conver UTM to lat lon
p2 = plot(lon,lat,'r');
scatter(lon(1), lat(1),'d','Cdata', [0,0,0])

[lat, lon] = utm2deg( xy2(:,1), xy2(:,2), repmat(utmzone,size(xy2,1),1) );  % Conver UTM to lat lon
s2 = scatter(lon, lat,'+','cdata',[1  0  1]);                                     % Plot Landmarks

hold off

% title('Locations of Side Scan Sonar Tacks')
% legend([p1,p2],{'Top','Bottom'})

clear p1 s1 p2 lat lon


%% Get feature Image Patches
% patches1 = GetPatches(image1, points1_filtered);
% patches2 = GetPatches(image2, points2_filtered);

% matches = ComparPatches(patches1, patches2);

%% FUNCTIONS ========================================================================================================================


%% Get feature Coordinates
function xy = GetFeatureCoords( points, depth, heading, midLine, utm)

xy = round( points );                                                       % X,Y coordinates of feature on SSS image [pixels]

theta = heading( xy(:,1) );                                                 % Heading of the vehicle coorisponding to the features [deg]

depth = depth(xy(:,1));                                                     % Height of the vehicle above sea floor [meters] 

dist = (xy(:,2) - midLine) .* 30/midLine;                                   % Range [meters]

dist = sqrt(dist.*dist - depth.*depth) .* sign(dist);                       % Horizontal distance of the feature [meters]

xy = utm( xy(:,1), : );                                                     % Location of the vehicle [utms "meters"]

theta = theta + 90;                                                         % Orientation of the feature [deg]

theta( theta < 0 )   = theta( theta < 0 ) + 360;                            % Keep  360 >= theta >= 0
theta( theta > 360 ) = theta( theta > 360 ) - 360;

xy(:,1) = xy(:,1) + dist .* sind( theta );                                  % Global position of the features
xy(:,2) = xy(:,2) + dist .* cosd( theta );

end

%% Compare Patches
function matches = ComparPatches(patches1, patches2)

matches = [];

for ii = 1: numel(patches1)
    
    im1 = patches1{ii};
    
    im1(isnan(im1)) = 0;
    
    for jj = 1: numel(patches2)
        im2 = patches2{jj};
        im2(isnan(im2)) = 0;
        
%         if numel(im1) > 1.5* numel(im2) || 1.5 * numel(im1) < numel(im2)
%             continue
%         end

        responce = conv2(im1, im2) ./numel(im2);
        
%         if numel(im1) < numel(im2)
%             responce = normxcorr2(im1, im2);
%         else
%             responce = normxcorr2(im2, im1);
%         end

%         if max(responce(:)) > 6500

        figure
        subplot(2,2,1)
        imshow(im1)
        
        subplot(2,2,2)
        imshow(im2)
        
        subplot(2,2,3)
        mesh(responce)
        w = sprintf("Responce: %f",max(responce(:)));
        title(w)
        
%     end
        
    end
end
end


%% Clear Features From Trench
function points = ClearTrench(points, midLine)

edge = 100;
ii = 1;

while ii <= points.Count
    
    if midLine + edge > points.Location(ii,2) && points.Location(ii,2) > midLine - edge
        points(ii,:) = [];
    else
        ii = ii+1;
    end
end
end


%% Merge Overlapping Features
function points = MergPoints(points)

ii = 1;
while ii <= points.Count
    
    jj = 1;
    while jj <= points.Count
        
        if jj == ii
            jj = jj+1;
            continue
        end
        
        x = points.Location(ii,1)- points.Location(jj,1);
        y = points.Location(ii,2)- points.Location(jj,2);
        dist = sqrt(x*x + y*y);
        
        if dist < 6 * (points.Scale(ii) + points.Scale(jj))
            
            points.Location(ii,1)  = mean([points.Location(ii,1),points.Location(jj,1)]);
            points.Location(ii,2)  = mean([points.Location(ii,2),points.Location(jj,2)]);
            points.Metric(ii)      = mean([points.Metric(ii),points.Metric(jj)]);
            points.Orientation(ii) = mean([points.Orientation(ii),points.Orientation(jj)]);
            points.Scale(ii)       = points.Scale(ii) + points.Scale(jj);
            points(jj,:) = [];
        else
            
        end
        
        jj = jj+1;
        
    end
    
    ii = ii+1;
    
end
end


%% Extract Image Patches
function patches = GetPatches(image, points, showPatches)

if nargin < 3
    showPatches = false;
end

patches{points.Count} = [];

for ii = 1: points.Count
    
    x = points.Location(ii,1);
    y = points.Location(ii,2);
    r = 6 * points.Scale(ii);
    
    rows = round( max(y-r,1): min(y+r, size(image,1)) );
    cols = round( max(x-r,1): min(x+r, size(image,2)) );
    
    patches{ii} = image(rows, cols);
    
    if  showPatches
        figure
        imshow(patches{ii})
    end
    
end
end


