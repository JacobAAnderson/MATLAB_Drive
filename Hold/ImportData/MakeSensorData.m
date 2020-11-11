close all
clear all
clc

% Create a timer object to queue the creation of data
%_Do not change these settings_________________________________________________________________
t = timer;
t.ExecutionMode = 'fixedRate';                  % Exicute TimerFcn at a fixed rate
t.TimerFcn = @MakeData;                         % Function to be exicuted by the timer
t.StopFcn = @(~,~) disp('Timer has Stopped');   % Tell the user that the timer has stopped
%______________________________________________________________________________________________


%__Change these settings_______________________________________________________________________
t.Period = 1;                                   % Execution Period in seconds
t.TasksToExecute = 10;                          % Number of exicutions before stopping
%______________________________________________________________________________________________


start(t);   % Start the timer

% TimerFcn -----------------------------------------------------------------------------------
function Data = MakeData(~,~)

% Generate imitation data
lon   = num2str(rand/100 - 107.8982160);
lat  =  num2str(rand/1000 + 37.2213890);
ph   = num2str(rand+7);
ORP  = num2str(rand*50+120);
DO   = num2str(rand*10);
temp = num2str(rand*10+20);
EC   = num2str(rand*800*1200);

% Concatenate the data
Data = [lon,',',lat,',',datestr(now,'yyyy-mm-dd HH:MM:SS'),';',ph,',',ORP,',',DO,',',temp,',',EC];

% Display data 
fprintf('Data: %s \n',Data)

end