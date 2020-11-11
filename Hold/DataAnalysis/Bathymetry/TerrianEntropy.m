close all
clear all
clc

ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(false);                                                     % Indicate whether new data should be selescted by the user

bathFile = ui.GetFile('bathymetry','mat');                                  % Get file path GUI
if isempty(bathFile), return, end                                           % End script if there isn't an input for it to use

bathy = LoadBathymetry(bathFile);                                      % Load Bathymetry Map


ti = Terrain_Info(bathy.AsMapp);
ti = ti.Eval_MapInfo(20);
ti.Map.InfoLayers
ti.Map.PlotLayer_Surf("Entropy");
ti.Map.PlotLayer_Surf("STD");
ti.Map.PlotLayer_Surf("Ave");

