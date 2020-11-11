% Get file Metadata in Unix system
clear all
close all
clc

% File of interest
fileName   = '/Users/jake/MATLAB-Drive/akaboticsGui/GeoReferancing.m';

% Set the Metadata search peramiters
lookingFor = 'FinderComment';       % Set this property in the MAC Finder / Get Info / Comment

% Call the OS to retrive Metadata fro the file
[status,cmdout] = system(['mdls ',fileName]);

% Check that the Metadata was returned succefully
if status
    disp('Metadat request failed')
    disp(status)
    return
end

% Pars out the data of interest
str = strsplit(cmdout,'\n');
str = str( contains(str,'FinderComment'));
str = strsplit(str{:},'"');
str = strsplit(str{2},';');
fileComment = str(1);
Requiernments = str(2);

% Display results
disp(fileComment)
disp(Requiernments)