% TAN POMDP

close all
clear all
clc

format long g

ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(false);                                                     % Indicate whether new data should be selescted by the user

bathFile = ui.GetFile('bathymetry','mat');                                  % Get file path GUI
if isempty(bathFile), return, end                                           % End script if there isn't an input for it to use

bathy = LoadBathymetry(bathFile);                                           % Load Bathymetry Map
clear bathFile ui                                                           % Clean up workspace

bathy = bathy.ReduceResolution(0.5);                                        % Reduce the Size of the Bathymetry Map

bathy = Terrain_Info(bathy.AsMapp);                                         % Evaluate Terrain Info
bathy = bathy.Eval_MapInfo(20);

bathy = bathy.Map;                                                          % Get just the map

bathy = bathy.Add_Layer('Policy');                                          % Add a policy layer to the map




start = [560943.53, 8092410.85];
goal  = [560775.86, 8094773.13];

% figure(fig)
% hold on 
% plot(start(1), start(2), '*r')
% plot(goal(1), goal(2), '*g')
% hold off


%%
clc
close all

fig = bathy.PlotLayer_Contour("STD");                                       % Show the map


idx = bathy.state2index(start);

X = bathy.index2state(idx)

neighbors = bathy.Neighbors(X, false)


figure(fig)
hold on 
plot(start(1), start(2), '*g')
plot(X(1), X(2), '*r')
plot(neighbors(:,1), neighbors(:,2), '*m')
hold off

