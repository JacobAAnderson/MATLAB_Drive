%% Demo surf2solid
close all
clear all
clc

load('C:\Users\jaanderson.FORTLEWIS\MATLAB Drive\_OutPutts\20171114_Bathymetry_LNH_Region_01.mat');
load('C:\Users\jaanderson.FORTLEWIS\MATLAB Drive\_OutPutts\10-Nov-2017_EcoMapperData.mat');

wqData = cat(1,data.wqData);
vehicle = cat(1,data.vehicle);

%%
fig = figure('Name','Plotting over Batymetry');
hold on
 surf(LON,LAT,Bathymetry,'EdgeColor','none'); 
 colormap(copper)
 caxis([-25 0])
 
 scatter3(vehicle(:,2),vehicle(:,1),vehicle(:,3),'.')
 
 hold off
