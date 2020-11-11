% Demo RRT
% Jacob Anderson
% 01/24/2020

close all
clear all
clc



%% create Environment
x = 0:0.01:20;
y = 0:0.01:10;

[X,Y] = meshgrid(x,y);

Z = 1./(1+exp(-(X-5)));


wall1 = [14, 10; 14, 5.75; 15, 5.75; 15,10];
wall2 = [14,  0; 14, 4.25; 15, 4.25; 15, 0];

map = Mapp(X, Y, Z, 'Environment');
map = map.Add_Obsticle(wall1);
map = map.Add_Obsticle(wall2);

map.PlotMap;
hold on

clear x X y Y Z wall1 wall2

start = [11.2673, 9.0816];
goal  = [17.8571, 1.4723];

plot(start(1), start(2), '*m')
plot(goal(1), goal(2), '*g')

rrt = RRT(map, 3);

path = rrt.Plan(start, goal, 1000, false);

plot(path(:,1), path(:,2), 'b')

hold off

