

close all
clear all
clc


ui = SavedUserInputs(mfilename);                                   % Instantiate the Saved User Inputs Class

ui = ui.NewData(true);                                                      % Indicate whether new data should be selescted by the user

file = ui.getFile('Data','mat');                                            % GUI to get the name and file path of a file
if isempty(file), return, end                                               % End script if there isn't an input for it to use

Bathymetry_Map = load(file);

figure(PlotBathymetry(Bathymetry_Map, 250, 50))                           % The last two inputs adjuxt the view angle
export_fig BathymetryMap.jpeg