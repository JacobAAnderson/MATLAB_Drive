


function rAd = Simulation1(num_trials, simulations, bathy, saveResults, wpFile1, wpFile2, count, plotIO)

[path, name, ext] = fileparts(saveResults);

filepath = fullfile(path, sprintf('%s_%d%s',name, count, ext));


dataTypes   = {'Err1', 'Err2', 'Time' 'Comms1', 'Coms2'};

rAd = RandomAssData([3000,num_trials], simulations, dataTypes);                 % Instanciate RandomAssData


% ___ Communications Planning _____________________________________________
auv1 = MR_AUV1( 'auv1', bathy, wpFile1, wpFile2);
auv2 = MR_AUV1( 'auv2', bathy, wpFile2, wpFile1);
dt = 0.5;

coms_policy1 = auv1.Planning.Plan(1600, dt, 1.25);
coms_policy2 = auv2.Planning.Plan(1600, dt, 1.25);

sz1 = numel(coms_policy1);
sz2 = numel(coms_policy2);

numcoms1 = sum(coms_policy1);
numcoms2 = sum(coms_policy2);

coms_policy1(end+1: 20000) = false;
coms_policy2(end+1: 20000) = false;



% ___ Run Simulation ______________________________________________________
for iter = 1 :num_trials
    for sim  = 1 : numel(simulations), fprintf("\n\nIteration: %i, Simulation: %s\n\n", iter, simulations{sim})
        
        % --- Fresh Robots ---
        clear auv1 auv2
        
        switch simulations{sim}
            case 'DR'
                auv1 = AUV_Basic( 'auv1', wpFile1);
                auv2 = AUV_Basic( 'auv2', wpFile2);
                
            otherwise
                auv1 = MR_AUV1( 'auv1', bathy, wpFile1, wpFile2);
                auv2 = MR_AUV1( 'auv2', bathy, wpFile2, wpFile1);
        end
        
        
        % --- Set up Random Ploicies ----
        switch simulations{sim}
            case 'rand05'
                rand_policy1 = false(20000,1);
                rand_policy2 = false(20000,1);
                a = 0.5;
                rand_policy1( randi(sz1, round(numcoms1*a), 1) ) = true;
                rand_policy2( randi(sz2, round(numcoms2*a), 1) ) = true;
                
            case 'rand1'
                rand_policy1 = false(20000,1);
                rand_policy2 = false(20000,1);
                a = 1;
                rand_policy1( randi(sz1, round(numcoms1*a), 1) ) = true;
                rand_policy2( randi(sz2, round(numcoms2*a), 1) ) = true;
                
            case 'rand2'
                rand_policy1 = false(20000,1);
                rand_policy2 = false(20000,1);
                a = 2;
                rand_policy1( randi(sz1, round(numcoms1*a), 1) ) = true;
                rand_policy2( randi(sz2, round(numcoms2*a), 1) ) = true;
                
            otherwise
        end
        
        
        
        
        % --- Reset Constants ---
        speed1 = 3;             % ASV1 Speed                        [m/s]
        speed2 = 3;             % ASV2 Speed                        [m/s]
        
        sim_time = 0;           % Simulation time                   [s]
        dt = 0.5;               % Simulation time step              [s]
        
        TIME(3000) = 0;         % Track Sim Time
        COMMS = false(3000,2); % Track Communications
        ii = 1;
        
        disp('Starting Simulation')
        tic
        
        % Run Simulation
        while toc < 1000
            
            % ____ Take Measurments _________________________________________
            auv1 = auv1.Add_SensorMeasurment("Depth",  -bathy.Depth( auv1.State.XYZ.GT(1:2)'), sim_time);
            auv1 = auv1.Add_SensorMeasurment("Compass", auv1.State.XYZ.GT(3), sim_time);
            
            auv2 = auv2.Add_SensorMeasurment("Depth",  -bathy.Depth( auv2.State.XYZ.GT(1:2)'), sim_time);
            auv2 = auv2.Add_SensorMeasurment("Compass", auv2.State.XYZ.GT(3), sim_time);
            
            
            % ___ Localization Update _____________________________________________
            dist1 = vecnorm(auv1.State.XYZ.GT(1:2)' - auv2.State.XYZ.GT(1:2)');
            
            switch simulations{sim}
                
                case 'DR' % Dear Reckoning Just Guesses
                    
                    
                case 'PF_EKF'
                    auv1 = auv1.UpdatePF(dt);
                    auv2 = auv2.UpdatePF(dt);
                    
                    COMMS(ii,:) = false(1,2);
                    
                    
                case 'PF_OneWay_Acoms'
                    auv1 = auv1.UpdatePF(dt);
                    
                    dist1 = vecnorm(auv1.State.XYZ.GT(1:2)' - auv2.State.XYZ.GT(1:2)');
                    acoms2 = {'acoms', dist1 + random('Normal', 0, 0.5), auv1.Nav.PF.X(1:2)', auv1.Nav.PF.Cov(1:2,1:2)};
                    auv2 = auv2.UpdatePF(dt, acoms2);
                    
                    COMMS(ii,:) = [true, false];
                    
                    
                    
                case 'PF_OneWay_Acoms_GT'
                    
                    auv1 = auv1.UpdatePF(dt);
                    
                    dist1 = vecnorm(auv1.State.XYZ.GT(1:2)' - auv2.State.XYZ.GT(1:2)');
                    acoms2 = {'acoms', dist1 + random('Normal', 0, 0.5), auv1.State.XYZ.GT(1:2)', [0.5, 0; 0, 0.5]};
                    auv2 = auv2.UpdatePF(dt, acoms2);
                    
                    COMMS(ii,:) = [true, false];
                    
                    
                    
                case 'PF_TwoWay_Acoms'
                    dist1 = vecnorm(auv1.State.XYZ.GT(1:2)' - auv2.State.XYZ.GT(1:2)');
                    
                    acoms1 = {'acoms', dist1 + random('Normal', 0, 0.5), auv2.Nav.PF.X(1:2)', auv2.Nav.PF.Cov(1:2,1:2)};
                    auv1 = auv1.UpdatePF(dt, acoms1);
                    
                    acoms2 = {'acoms', dist1 + random('Normal', 0, 0.5), auv1.Nav.PF.X(1:2)', auv1.Nav.PF.Cov(1:2,1:2)};
                    auv2 = auv2.UpdatePF(dt, acoms2);
                    
                    COMMS(ii,:) = true(1,2);
                    
                 
                    
                case 'PF_ConstraintOnly_Acoms'
                    % ___ Communications Planning _________________________________________
                    [auv1.Planning, comms1] = auv1.Planning.Update_Belief(dt, auv1.Nav.PF);
                    [auv2.Planning, comms2] = auv2.Planning.Update_Belief(dt, auv2.Nav.PF);
                    
                    COMMS(ii,:) = [comms1, comms2];
                    
                    % ___ A-Comms _________________________________________________________
                    dist1 = vecnorm(auv1.State.XYZ.GT(1:2)' - auv2.State.XYZ.GT(1:2)');
                    
                    if comms2
                        acoms1 = {'acoms', dist1 + random('Normal', 0, 0.5), auv2.Nav.PF.X(1:2)', auv2.Nav.PF.Cov(1:2,1:2)};
                        auv1 = auv1.UpdatePF(dt, acoms1);
                    else
                        auv1 = auv1.UpdatePF(dt);
                    end
                    
                    
                    if comms1
                        acoms2 = {'acoms', dist1 + random('Normal', 0, 0.5), auv1.Nav.PF.X(1:2)', auv1.Nav.PF.Cov(1:2,1:2)};
                        auv2 = auv2.UpdatePF(dt, acoms2);
                    else
                        auv2 = auv2.UpdatePF(dt);
                    end
                    
                 
                    
                case 'Policy'
                    if coms_policy2(ii)
                        acoms1 = {'acoms', dist1 + random('Normal', 0, 0.5), auv2.Nav.PF.X(1:2)', auv2.Nav.PF.Cov(1:2,1:2)};
                        auv1 = auv1.UpdatePF(dt, acoms1);
                    else
                        auv1 = auv1.UpdatePF(dt);
                    end
                    
                    
                    if coms_policy1(ii)
                        acoms2 = {'acoms', dist1 + random('Normal', 0, 0.5), auv1.Nav.PF.X(1:2)', auv1.Nav.PF.Cov(1:2,1:2)};
                        auv2 = auv2.UpdatePF(dt, acoms2);
                    else
                        auv2 = auv2.UpdatePF(dt);
                    end
                    
                    COMMS(ii,:) = [coms_policy1(ii), coms_policy2(ii)];
                    
                    
                    
                    
                case {'rand05', 'rand1', 'rand2'}
                    dist1 = vecnorm(auv1.State.XYZ.GT(1:2)' - auv2.State.XYZ.GT(1:2)');
                    
                    if rand_policy2(ii)
                        acoms1 = {'acoms', dist1 + random('Normal', 0, 0.5), auv2.Nav.PF.X(1:2)', auv2.Nav.PF.Cov(1:2,1:2)};
                        auv1 = auv1.UpdatePF(dt, acoms1);
                    else
                        auv1 = auv1.UpdatePF(dt);
                    end
                    
                    
                    if rand_policy1(ii)
                        acoms2 = {'acoms', dist1 + random('Normal', 0, 0.5), auv1.Nav.PF.X(1:2)', auv1.Nav.PF.Cov(1:2,1:2)};
                        auv2 = auv2.UpdatePF(dt, acoms2);
                    else
                        auv2 = auv2.UpdatePF(dt);
                    end
                    
                    COMMS(ii,:) = [rand_policy1(ii), rand_policy2(ii)];
                    
                    
                otherwise
                    warning('Simulation Not Found')
            end
            
            
            
            % ___ ASV Operations __________________________________________________
            auv1 = auv1.Navigate(speed1, dt , 5);
            auv2 = auv2.Navigate(speed2, dt , 5);
            
            
            
            % ___ Time Keeping ____________________________________________
            TIME(ii) = sim_time;
            ii = ii+1;
            sim_time = sim_time +dt;
            
            
            
            % ____ Check for goal _________________________________________________
            if any( [auv1.Goal_Achived(0.5), auv1.Nav.Done, auv2.Goal_Achived(0.5), auv2.Nav.Done])
                
                % Run Stats ===========================================================
                stats1 = auv1.PathStates(false);
                stats2 = auv2.PathStates(false);
                
                % Save Data
                for jj = 1 : ii-1
                    rAd = rAd.Add_Result(iter, simulations{sim}, stats1.err(jj), stats2.err(jj), TIME(jj), COMMS(jj,1), COMMS(jj,2)); % Add data to the respective set
                end
                
                
                save(filepath, 'rAd')
                
                
                % =====================================================================
                
                break
            end
            
        end
        
        
        
        fprintf('\nSimulation Time: %f [min]\nRun Time: %f [min]\n', sim_time/60, toc/60)
        disp('_____________________________________________________')
        
        rAd = rAd.ResetIndex;
        
        clear dist1 comms1 comms2 acoms1 acoms2 dt speed1 speed2 sim_time
        
        
        if plotIO
            fig = bathy.Plot_3DModel(0, 90);                                % Show AUV tracks on the Bathymetry Map
            
            auv1.Nav.PlotWaypoints(fig);
            auv2.Nav.PlotWaypoints(fig);
            
            auv1.PlotPaths(fig);
            auv2.PlotPaths(fig);
            drawnow
        end
        
        
        
    end
end


% Joint AUV Error
flds = fields(rAd.Data);

for ii = 1: numel(flds)
    err = rAd.Data.(flds{ii}).Err1 + rAd.Data.(flds{ii}).Err2;
    rAd.Data.(flds{ii}).JointErr = err;
end

rAd = rAd.Eval_Stats;
save(filepath, 'rAd')


end
%___________________________________________________________________________________________________






