%% Matlab file to load data from a ROS bagfile created by a heron USV
close all
clear all
clc

%%
tic
[fileName, filePath] = uigetfile('*.bag','Select .bag file','C:\Users\jaanderson.FORTLEWIS\_WorkingData');
if fileName == 0 
    return
end

file = fullfile(filePath,fileName);
bag = rosbag(file);
toc

fprintf('\n\nDuration of Mission: %f Seonds? \n\n',(bag.EndTime-bag.StartTime))


%% Get Latitude and Longitude
nav_fix  = select(bag, 'Time',[bag.StartTime bag.StartTime+60], 'Topic', '/navsat/fix');    % Select the messages to read from the .bag file
nav_fix_msg = readMessages(nav_fix);
ts_nav_fix  = timeseries(nav_fix);                                                             % Creat a time series from the messages
lat_lon = ts_nav_fix.Data(:,1:2);                                                              % Extract latatude on longitude



%% Get GPS NMEA sentance
nav_nmes     = select(bag, 'Time',[bag.StartTime bag.StartTime+60], 'Topic', '/navsat/rx');
msg_nav_nmea = readMessages(nav_nmes);
nmea = cellfun(@(A) A.Sentence_,msg_nav_nmea,'UniformOutput',false)


%% Other Nav Data
% nav_stat = select(bag, 'Time',[bag.StartTime bag.StartTime+60], 'Topic', '/navsat_status');
% nav_time = select(bag, 'Time',[bag.StartTime bag.StartTime+60], 'Topic', '/navsat/time_reference');
% 
% msg_nav_stat = readMessages(nav_stat);
% msg_nav_time = readMessages(nav_time);
% 
% ts_nav_Stat = timeseries(nav_stat);
% ts_nav_Time = timeseries(nav_time);


%% Get Sonde Data
% rosgenmsg('C:\Users\jaanderson.FORTLEWIS\Documents\HeronCode\ROS_PackagesForHeron');
% sonde = select(bag, 'Time',[bag.StartTime bag.StartTime+100], 'Topic', '/sonde');
% sondeData = readMessages(sonde); 
% ,'Format', 'MM.dd.yyyy HH:mm:ss.SS');

