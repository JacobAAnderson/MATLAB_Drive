% Ocean Color Data From NASA =====================================================================================================
close all
clear all
clc


%% Build URL for Ocean Data -----------------------------------------------------------
url = 'https://tds.hycom.org/thredds/dodsC/GLBy0.08/latest';

finfo = Extract_NetCF_File(url);