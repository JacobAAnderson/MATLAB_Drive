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

clear x X y Y Z wall1 wall2

x = [14.5, 8];
sig = [9,6;6,9];

tf = map.Collition(x);

[hit, dev] = map.CollitionDist(x, sig);



