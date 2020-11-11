%% Get a geo-referanced mapp from an online Web Map Server
close all
clear all
clc


%% Bounding Box: Lower Left Corner and Uper Rigth Corner
lon = [-107.8929390347302,  -107.9508463873537];    % Min and Max longitudes
lat = [  37.20918907757923,   37.24578625143064];   % Min and Max latitudes



%% Get Geottiff
[geoImage, geoData] = GetMapfromWMS( lon, lat );



%% Disply the Map 
geoshow(geoImage, geoData)

