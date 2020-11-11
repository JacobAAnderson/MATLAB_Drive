% Demo Geotiff Class
% Jake Anderson
% 9/25/2019

close all
clear all
clc


ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(true);                                                      % Indicate whether new data should be selescted by the user


file = ui.GetFile('geotiff','tiff');                                        % GUI to get the name and file path of a file
if isempty(file), return, end                                               % End script if there isn't an input for it to use

gt = Geotiff;

gt = gt.FromFile(file);
figure
gt.Show;


lon = [-107.9113366574420,  -107.9073118934940];                            % Min and Max longitudes
lat = [  37.2387203649409,    37.2387286862432];                            % Min and Max latitudes
zoom = 18;

gt = gt.FromGoogle(mean(lat), mean(lon), zoom);
figure
gt.Show;



gt = gt.FromWMS(lat, lon);
figure
gt.Show;
