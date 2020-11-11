
close all
clear all
clc

ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(false);                                                     % Indicate whether new data should be selescted by the user

bathFile = ui.GetFile('bathymetry','mat');                                  % Get file path GUI
if isempty(bathFile), return, end                                           % End script if there isn't an input for it to use

bathymetry = LoadBathymetry(bathFile);                                      % Load Bathymetry Map


wp = Mapp(bathymetry.Easting, bathymetry.Northing, bathymetry.Elevation);


X = [362000, 3701500];

bathymetry.Depth(X)
wp.Measure(X)