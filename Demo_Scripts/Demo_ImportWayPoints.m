%% Demo scrip to call ImportWaypoints function
close all
clear all
clc

wayPoint_filePath = 'C:\Users\jaanderson.FORTLEWIS\Desktop\';
wayPoint_fileName1 = 'waypoints.txt';
%wayPoint_fileName2 = 'waypoints2.txt';

wayPoints1 = ImportWayPoints(wayPoint_filePath,wayPoint_fileName1,{' '});
%wayPoints2 = ImportWayPoints(wayPoint_filePath,wayPoint_fileName2,{' '});

% wayPoints1 = flipud(wayPoints1); % Flips the waypoint linst upside down
wayPoints1 = fliplr(wayPoints1); % Switch order of Lat and Long


%% Example outputs
%

fid=fopen('WayPoints.txt','a');
fprintf(fid,'%3.12f,%2.13f\r\n', wayPoints1');
fclose(fid);
