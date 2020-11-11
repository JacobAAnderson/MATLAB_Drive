close all
clear all
clc
format long g

ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(false);                                                     % Indicate whether new data should be selescted by the user

bathFile = ui.GetFile('bathymetry','mat');                                  % Get file path GUI
if isempty(bathFile), return, end                                           % End script if there isn't an input for it to use

bathymetry = load(bathFile);                                                % load Bathymetry map
bathymetry = bathymetry.bathy_map;                                          % Get rid of loading structue
clear bathFile


wf = WaterFeature3(bathymetry);
fig_map = wf.Map.PlotMap_Contour;

% wf = wf.Find_Gradiant;
% fig_grad = wf.PlotGrad;
% 
% wf.Plot_MaxMin([0,90], fig_map);
% 
% wf = wf.Fit_Model;
% wf.PlotModel;



%% Move Stuff Around

clc
wf = wf.SetUncertainty('Speed',     'Normal', -0.01, 0.005);
wf = wf.SetUncertainty('Compass',   'Normal', -5,   10);
wf = wf.SetUncertainty('Spreading', 'Normal',  0,    0.2);

% Move Feature
loops = 30;
movie = Fig2Movie('Water Feature Sim', loops);
for ii = 1: loops
    
    speed = 0.2;
    theta = 320;
    dt = 1;
    
    wf = wf.Move(speed, theta, dt);
    wf = wf.MakeMap;
    wf.PlotMap([0,90],fig_map);
    
    movie = movie.Add_Frame(fig_map);
    
end


% Save file GUI ---------------------------------------------------------------------------------------------------------
ui = ui.NewData(true);

saveFile = ui.SaveFile('water_Feature','mp4');                              % GUI to get file path to save variable / object
if isempty(saveFile), return, end                                           % End script if there isn't an input for it to use

%%
movie.MakeVideo(saveFile, 5)


disp("Done!!")

