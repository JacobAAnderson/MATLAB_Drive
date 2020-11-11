close all
clear all
clc

ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(true);                                                     % Indicate whether new data should be selescted by the user

% Get Bathymetry Map
bathFile = ui.GetFile('bathymetry','mat');                                  % Get file path GUI
if isempty(bathFile), return, end                                           % End script if there isn't an input for it to use

bathy = LoadBathymetry(bathFile);                                           % Load Bathymetry Map

bathy = bathy.ReduceResolution(0.5);

% Make terrain info maps
ti = Terrain_Info(bathy.AsMapp);
ti = ti.Eval_MapInfo(20);


%% Access info
ti.Map.InfoLayers
ti.Map.Measure(bathy.Center, "Entropy")
ti.Map.Measure(bathy.Center, "STD")
ti.Map.Measure(bathy.Center, "Ave")


% Display the information maps
ti.Map.PlotLayer_Surf("Entropy");
ti.Map.PlotLayer_Surf("STD");
ti.Map.PlotLayer_Surf("Ave");


