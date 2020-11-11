% Template MATLAB code for reading data from a private channel, analyzing
clear all
close all
clc



readChannelID = 411364;

readAPIKey = 'LAVZF4O5ZJ6NTC4K';

% TODO - Replace the [] with channel ID to write data to:
%writeChannelID = [];
% TODO - Enter the Write API Key between the '' below:
%writeAPIKey = '';

%% Read Data %%
dataTable = thingSpeakRead(readChannelID, 'ReadKey', readAPIKey,'Fields',1,'OutputFormat','table','NumMinutes',10)


%% Analyze Data %%

data = dataTable{:,2}               % Extract Data from the Table

splitData = strsplit(data{1},';')   % Split the CSV

% Data Header:
% "UTC Date [ddmmyy], UTC Time [hhmmss.sss], Lat [ddmm.mmmm], Lon [ddmm.mmmm];Battery [V]; Array1: pH, DO [mg/L], ORP [mV], EC [uS/cm], TDS [ppm], Sal, SG, temp [C], pres [bar], depth [m], temp [c]; Array2: pH, DO [mg/L], ORP [mV], EC [uS/cm], TDS [ppm], Sal, SG, temp [C], pres [bar], depth [m], temp [c] ";

%% Write Data %%
%thingSpeakWrite(writeChannelID, analyzedData, 'WriteKey', writeAPIKey);

