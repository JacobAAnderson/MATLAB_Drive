%% Filter Waypoints
%  Isolates a set of waypoint inside a specified region
close all
clear all
clc


%% Get working data
% Use Geotiff for visual feedback
disp('Choose Geotiff for data visualization')
disp('  This is opptioanl')
[geotiffFileName,geotiffFilePath] = uigetfile('*.tif*','Choose Geotiff');

if geotiffFileName ~= 0
    geoTifff = fullfile(geotiffFilePath,geotiffFileName);
    [geoImage,geoData] = geotiffread(geoTifff);
    figure('Name','Map')
    map = geoshow(geoImage,geoData);
    hold on
end

% Get waypoint files
disp('Select waypoints to filter')
disp('  i.e: The waypoints that you want a subset of')
[innerFileName,innerFilePath] = uigetfile('*.txt','Select Waypoint to filter');


disp('Select wapoints of the region of interest')
disp('  i.e: The waypoint that define the area of interest')
[outerFileName,outerFilePath] = uigetfile('*.txt','Select Waypoint to filter',innerFilePath);


%% Read in waypoint file
wayPoints  = ImportWayPoints(innerFilePath,innerFileName,'\n');
edgePoints = ImportWayPoints(outerFilePath,outerFileName,'\n');

% Make sure waypoints are Lon - Lat
if wayPoints(1,1) > 0
    wayPoints = fliplr(wayPoints);
end

if edgePoints(1,1) > 0
    edgePoints = fliplr(edgePoints);
end


% Plot waypoint files for visual feedback
p1 = plot(wayPoints(:,1),wayPoints(:,2));
legend(p1,'Waypoints','Location','north')
pause(2)

p2 = plot(edgePoints(:,1),edgePoints(:,2));
legend([p1,p2],{'Waypoints','Enclosing Area'},'Location','north')
pause(2)


% Filter the waypoints
in = inpolygon(wayPoints(:,1),wayPoints(:,2),edgePoints(:,1),edgePoints(:,2));
wayPoints = wayPoints(in,:);

% Dispaly the new waypoints
delete(p1)

p1 = plot(wayPoints(:,1),wayPoints(:,2));
legend([p1,p2],{'New Waypoints','Enclosing Area'},'Location','north')
pause(2)


% Save the new waypoint list as a txt file --> lat lon
[saveFileName,savePathName] = uiputfile('.txt','Save new waypoints');
savefile = fullfile(savePathName, saveFileName);


%%
fid = fopen(savefile,'w');
fprintf(fid,'%3.13f,%4.13f \n',[wayPoints(:,2),wayPoints(:,1)]');
fclose(fid);


