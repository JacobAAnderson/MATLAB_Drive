%%% Import Sensor Buoy Data
clear all
close all
clc

%% Get Data
% Get File info ------------------------------------------------------------------------------
[data_name, data_filepath] = uigetfile('*.csv*','Select Data');

if data_name == 0
    return
end

% Read data file -----------------------------------------------------------------------------
dataFile = fullfile(data_filepath, data_name);

delimiter = '"';

fileID = fopen(dataFile,'r');
dataArray = textscan(fileID,'%s'); %, 'Delimiter', delimiter, 'HeaderLines' ,1, 'ReturnOnError', false);
fclose(fileID);

dataArray = dataArray{1};

header = strsplit(dataArray{1},',');

% Pre-alocate memory to data variables -------------------------------------------------------------------------------------
num = length(dataArray);

time(num,1)      = datetime('20140304','InputFormat','yyyyMMdd');
longitude(num,1) = 0;
latitude(num,1)  = 0;
pH(num,2)        = 0;
DO(num,2)        = 0;
ORP(num,2)       = 0;
EC(num,2)        = 0;
Sal(num,2)       = 0;
SG(num,2)        = 0;
pressure(num,2)  = 0;
depth(num,2)     = 0;
temp(num,2)      = 0;

index = 1;
lineErros = 0;
wrongDate = 0;
for ii = 2 : num
    try
        line = strsplit(dataArray{ii},'"');
        
        Array1 = strsplit(line{2},',');
        Array2 = strsplit(line{4},',');
        GPS    = strsplit(line{6},',');
        timestamp = strsplit(line{1},{'T','Z'});
        
%        time(index) =  datetime(strcat(GPS{1},GPS{2}),'InputFormat','ddMMyyHHmmss.SSS');        % GPS Time Stamp
        time(index) = datetime([timestamp{1},' ',timestamp{2}],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');       % Graylog time stamp
        
        if time(index) < datetime('20180420','InputFormat','yyyyMMdd')
            wrongDate = wrongDate +1;
            continue
        end
        
        longitude(index,1) = str2double(GPS{3});
        latitude(index,1)  = str2double(GPS{4});
        pH(index,:)        = [ str2double(Array1{ 1}), str2double(Array2{ 1}) ];
        DO(index,:)        = [ str2double(Array1{ 2}), str2double(Array2{ 2}) ];
        ORP(index,:)       = [ str2double(Array1{ 3}), str2double(Array2{ 3}) ];
        EC(index,:)        = [ str2double(Array1{ 5}), str2double(Array2{ 5}) ];
        Sal(index,:)       = [ str2double(Array1{ 6}), str2double(Array2{ 6}) ];
        SG(index,:)        = [ str2double(Array1{ 8}), str2double(Array2{ 8}) ];
        pressure(index,:)  = [ str2double(Array1{ 9}), str2double(Array2{ 9}) ];
        depth(index,:)     = [ str2double(Array1{10}), str2double(Array2{10}) ];
        temp(index,:)      = [ str2double(Array1{11}), str2double(Array2{11}) ];
        
        index = index+1;
        
    catch
        lineErros = lineErros + 1;
    end
end

% Trunckate unused end of array
time(index:num,:)      = [];
longitude(index:num,:) = [];
latitude(index:num,:)  = [];
pH(index:num,:)        = [];
DO(index:num,:)        = [];
ORP(index:num,:)       = [];
EC(index:num,:)        = [];
Sal(index:num,:)       = [];
SG(index:num,:)        = [];
pressure(index:num,:)  = [];
depth(index:num,:)     = [];
temp(index:num,:)      = [];

% Sort rows by date and time
[time, index] = sortrows( time );

longitude = longitude(index,:);
latitude  = latitude(index,:);
pH        = pH(index,:);
DO        = DO(index,:);
ORP       = ORP(index,:);
EC        = EC(index,:);
Sal       = Sal(index,:);
SG        = SG(index,:);
pressure  = pressure(index,:);
depth     = depth(index,:);
temp      = temp(index,:);

fprintf("Line Errors: %d \n Date Errors: %d", lineErros, wrongDate);



%% Do stuff with data

figure('Name','pH','NumberTitle','off')
hold on
plot(time,pH(:,1));
plot(time,pH(:,2));
hold off
title('pH at Depth')
legend({'Array1', 'Array2'},'Location','northeastoutside')
xlabel('Date')
ylabel('pH')
set(get(gca,'ylabel'),'rotation',0,'HorizontalAlignment','right')


figure('Name','DO','NumberTitle','off')
hold on
plot(time,DO(:,1));
plot(time,DO(:,2));
hold off
title('Dissolved Oxygen at Depth [mg/L]')
legend({'Array1', 'Array2'},'Location','northeastoutside')
xlabel('Date')
ylabel('mg/L')
set(get(gca,'ylabel'),'rotation',0,'HorizontalAlignment','right')


figure('Name','ORP','NumberTitle','off')
hold on
plot(time,ORP(:,1));
plot(time,ORP(:,2));
hold off
title('ORP at Depth in [mV]')
legend({'Array1', 'Array2'},'Location','northeastoutside')
xlabel('Date')
ylabel('mV')
set(get(gca,'ylabel'),'rotation',0,'HorizontalAlignment','right')


figure('Name','EC','NumberTitle','off')
hold on
plot(time,EC(:,1));
plot(time,EC(:,2));
hold off
title('Conductivity at Depth in [\muS/cm]')
legend({'Array1', 'Array2'},'Location','northeastoutside')
xlabel('Date')
ylabel('\muS/cm')
set(get(gca,'ylabel'),'rotation',0,'HorizontalAlignment','right')


figure('Name','temp','NumberTitle','off')
hold on
plot(time,temp(:,1));
plot(time,temp(:,2));
hold off
title('Temperature at Depth [\circC]')
legend({'Array1', 'Array2'},'Location','northeastoutside')
xlabel('Date')
ylabel('\circC')
set(get(gca,'ylabel'),'rotation',0,'HorizontalAlignment','right')


figure('Name','depth','NumberTitle','off')
hold on
plot(time,depth(:,1));
plot(time,depth(:,2));
hold off
title('Sensor Node Depths [Meters]')
legend({'Array1', 'Array2'},'Location','northeastoutside')
xlabel('Date')
ylabel('Meters')




