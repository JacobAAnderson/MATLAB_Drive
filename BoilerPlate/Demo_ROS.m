%% Demo ROS Script
%  Demo Script on using the Robot Operating System with matlab
%  Use with Turtlebot Stage Simulation
%  roslaunch turtlebot_stage turtlebot_in_stage.launch

close all
clear all
clc

%% Start up
% Connect to an external ros master
hostName = 'Me';                                                            % Name of the computer hosting ROS Master
rosinit(hostName);                                                          % Conect to ros master

% Start a ros master
% rosinit

%% Initialize Publishers and subscribers
chatterpub = rospublisher('/chatter', 'std_msgs/String');                   % Publish Data
pause(2)                                                                    % Wait to make sure publisher is registered

laser  = rossubscriber('/scan');                                            % Subscribe to a topic
listen = rossubscriber('/chatter', @listenrFunction);                       % Subscribe to a topic with a callback function


%% Recive Data
scandata = receive(laser,10);                                               % Get mesage

figure('Name','Laser Data')                                                 % Plot Data
plot(scandata,'MaximumRange',7)


%% Publis Data
disp('Getting Ready to publish data to "/chatter" topic')
disp('  rostopic echo chatter')
disp('Paused, Press and button to continue')
pause

chattermsg = rosmessage(chatterpub);                                        % Create a message
for ii = 0: 10
    
chattermsg.Data = sprintf('hello world %d', ii);                                            % Assigne data
send(chatterpub,chattermsg)                                                 % Send data
pause(1)                                                                    % Allow message to be transmited

end

%% Shutdown ROS
rosshutdown


%% Call Back Function
function listenrFunction(subscriber, msgs)

fprintf("\n\n\n\n____ New Message _____________________________________________________________\n")
disp(subscriber(:))
disp(msgs(:))

end
