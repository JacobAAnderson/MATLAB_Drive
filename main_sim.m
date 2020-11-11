% Robot Simulator

close all
clear all
clc

% ___ Load Data ___________________________________________________________
ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(false);                                                     % Indicate whether new data should be selescted by the user

bathFile = ui.GetFile('bathymetry','mat');                                  % Get file path to Bathymetry Map
if isempty(bathFile), return, end                                           % End script if there isn't an input for it to use

wpFile = ui.GetFile('waypoints','txt');                                     % Get file path to waypoint list
if isempty(wpFile), return, end                                             % End script if there isn't an input for it to use


saveResults = ui.SaveFile('Results','mat');                                 % Get File Path to save Results
if isempty(saveResults), return, end                                        % End script if there isn't an input for it to use


saveStats = ui.SaveFile('Stats','tex');                                     % Get File Path to save statistics from results
if isempty(saveStats), return, end                                          % End script if there isn't an input for it to use


% ui = ui.NewData(true);
% saveFile = ui.SaveFile('Sim_Run','mp4');                                    % Get File Path to save movie
% if isempty(saveFile), return, end                                           % End script if there isn't an input for it to use


bathy = LoadBathymetry(bathFile);                                           % Load Bathymetry Map
clear bathFile      

trials = 20;
rad = RandomAssData(trials, ['env1'; 'env2'; 'env3'], ['Erro'; 'Perc'; 'Time']);

for iter = 1 : trials
    
    fprintf("Iteration: %i\n", iter)
    
    for env = 1:3
        
        switch env
            case 1, wf = WaterFeature1(bathy, false);                                           % Create a water Feature with the bathymetry mas as the foot print
            case 2, wf = WaterFeature2(bathy, false);
            case 3, wf = WaterFeature3(bathy, false);
            otherwise
                continue
        end
        % ___ Create ASV __________________________
        asv = ASV1(bathy);                                                          % Instanciate an ASV with the bathymetry maps global referance fram
        asv = asv.LoadWaypoints(wpFile);                                            % Give the ASV a set of waypoints to follow
        
        % ___ Create AUV __________________________
        auv = AUV1(bathy);
        
        
        % ___ Diplay Initial Setup ________________
        % fig = bathymetry.Plot_3DModel(-97, 28);                                   % Plot Bathymetry Model
        % fig = wf.Map.PlotMap_Contour; %(fig, 0.7);
        % asv.PlotLocation(fig);
        % asv.Nav.PlotWaypoints(fig);
        % pause(2)
        
        
        %% Run Simulation
        % movie = Fig2Movie('Robot Sim Recorder');                                    % Create Movie Recorder
        
%         run_wf = false;
%         wf_speed = 0;                                                               % Water Feature Speed               [m/s]
%         wf_theta = 140;                                                             % Water Feature Direction of Travel [deg]
        
        run_asv = true;
        asv_r = 3;                                                                  % ASV Waypoint Achievment Radius    [m]
        asv_speed = 3;                                                              % ASV Speed                         [m/s]
        
        run_auv = false;
        auv_r = 3;                                                                  % ASV Waypoint Achievment Radius    [m]
        auv_speed = 3;                                                              % ASV Speed                         [m/s]
        
        
        sim_time = 0;                                                               % Simulation Time                   [s]
        dt = 0.5;                                                                   % Simulation Time Step              [s]
        
        tic
        while toc < 1000
            
            % ___ Move the Water Feature __________________________________________
            %     if run_wf
            %         wf = wf.Move(wf_speed, wf_theta, dt);
            %         wf = wf.MakeMap;
            %         fig = wf.Map.PlotMap(fig, 0.7);
            %     end
            
            
            %___ ASV Operations __________________________________________________
            if run_asv
                asv = asv.Navigate(asv_speed, dt , 3);
                %         asv = asv.Add_SensorMeasurment("Bathy", bathy.Depth( asv.State.XYZ.GT(1:2)), sim_time);
                asv = asv.Add_SensorMeasurment("Salt",  wf.Map.Measure( asv.State.XYZ.GT(1:2)),   sim_time);
                %         asv = asv.Step(sim_time);
                %         asv.PlotLocation(fig);
                %         asv.PlotPaths(fig)
                
            end
            
            
            if asv.Nav.Done
                
                %         asv.PlotPaths(fig)
                %         movie = movie.Add_StillFrame(fig, 2);
                
                run_auv = true;
                run_asv = false;
                
                asv.Nav.Done = false;
                %         asv.PlotPaths(fig)
                asv = asv.PredictFeatures("Salt", false);
                
                %         fig2 = asv.Models.Salt.PlotMap;
                %         movie = movie.Add_StillFrame(fig2, 2);
                %
                %         fig2 = asv.Models.Salt.Plot_MaxMin([0,90], fig2);
                %         movie = movie.Add_StillFrame(fig2, 2);
                
                
                auv = auv.ReciveWaterPropertyMap("Salt", asv.SendWaterPropertyMap("Salt"));
                
                
                %         fig = auv.Models.Salt.Map.PlotMap_Contour;
                
            end
            
            
            % ___ AUV Operations __________________________________________________
            if run_auv
                auv = auv.Navigate(auv_speed, dt , 3);
                auv = auv.Add_SensorMeasurment("Bathy", bathy.Depth( auv.State.XYZ.GT(1:2)), sim_time);
                auv = auv.Add_SensorMeasurment("Salt",  wf.Map.Measure( auv.State.XYZ.GT(1:2)),   sim_time);
                %         auv = auv.Step(sim_time);
                %         auv.PlotLocation(fig);
                %         auv.PlotPaths(fig)
                
                %         movie = movie.Add_Frame(fig);                                           % Record Figure
            end
            
            
            sim_time = sim_time +dt;                                                % time keeping
            
            if auv.Nav.Done, break; end
            
        end
        
        hold off
        
        
        %% Wrap up
        % auv.PlotPaths(fig)
        % auv.Models.Salt.Plot3DModel(asv.Map)
        
        auv = auv.PredictFeatures("Salt", false);
        
%         fig = auv.Models.Salt.Map.PlotMap_Contour;                                  % Add Results to the film
        % movie = movie.Add_StillFrame(fig, 2);
        %
        % fig = wf.Map.PlotMap_Contour;
        % movie = movie.Add_StillFrame(fig, 2);
        %
        % FrameRate = 30;                                                             % Frames per second
        % movie.MakeVideo(saveFile, FrameRate)                                        % Make Movie
        
        
        %% Calculate RMSE
        A = auv.Models.Salt.Map.Map;
        B = wf.Map.Map;
        
        B = imresize(B, [size(A)]);
        
        mse = MSError(A, B);
        
        A(~isnan(A)) = 0;
        nom = MSError(A, B);
        
        per =10 * mse/nom;
        
        fprintf("MRSE: %f, per: %f\n", mse, per)
        
        
        
        switch env
            case 1, rad = rad.Add_Result('env1', mse, per, sim_time);                                          % Create a water Feature with the bathymetry mas as the foot print
            case 2, rad = rad.Add_Result('env2', mse, per, sim_time);
            case 3, rad = rad.Add_Result('env3', mse, per, sim_time);
            otherwise
                continue
        end
        
        clear asv auv mes nom A B sim_time
        
    end
end

save(saveResults, 'rad')


%%
rad = rad.Eval_Stats;

rad.Write_Stats(saveStats, 'Erro', '%1.5f')
rad.Write_Stats(saveStats, 'Perc', '%1.4f')




