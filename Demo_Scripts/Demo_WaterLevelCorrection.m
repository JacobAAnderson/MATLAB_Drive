% Script to call WaterLevelCorrection function

close all
clear all
clc
format longG

filePath = 'C:\Users\jaanderson.FORTLEWIS\MATLAB Drive\_WorkingData\';
txtName = 'LNH_WaterLevelOffset.txt';

[bathy,DATA] = ImportWaterLevel(filePath,txtName);