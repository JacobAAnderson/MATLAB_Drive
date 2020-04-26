
% simulations = {'Policy', 'DR', 'PF_EKF', 'rand05', 'rand1', 'rand2'};


function rAd = Simulation2(num_trials, num_auv, simulations, bathy, saveResults, wpFile, count, plotIO)

[path, name, ext] = fileparts(saveResults);

filepath = fullfile(path, sprintf('%s_%d%s',name, count, ext));

for vv = num_auv:-1:1
    dataTypes{vv} = sprintf('Err%d',vv);
    dataTypes{num_auv+1+vv} = sprintf('Comms%d',vv);
end
dataTypes{num_auv+1} = 'Time';

rAd = RandomAssData([3000,num_trials], simulations, dataTypes);                 % Instanciate RandomAssData


% ___ Communications Planning _____________________________________________
for vv = num_auv:-1:1
   
   in = true(1,num_auv);
   in(vv) = false;
   auv(1:vv) = MR_AUV1( sprintf('auv%d',vv), bathy, wpFile{vv}, wpFile{in});
   
end


if plotIO
    fig = bathy.Plot_3DModel(0, 90);                                % Show AUV tracks on the Bathymetry Map
    for jj = 1:num_auv, auv(jj).Nav.PlotWaypoints(fig); end
    drawnow
end



coms_policy = zeros(num_auv, 2000);

fprintf('\n\nStarting Communications Planning\n\n')

dt = 0.5;
for vv = num_auv: -1: 1
    
    policy = auv(vv).Planning.Plan1(1600, dt, 1.25);
    coms_policy(vv, 1:size(policy,2) ) = policy;
 
    sz(vv) = numel(policy);
    numcoms(vv) = sum(policy);

end
fprintf('\n\nFinished Communications Planning\n\n')



% ___ Run Simulation ______________________________________________________
for iter = 1 :num_trials
    for sim  = 1 : numel(simulations), fprintf("\n\nIteration: %i, Simulation: %s\n\n", iter, simulations{sim})
        
        % --- Fresh Robots ---
        clear auv
                
        for vv = num_auv:-1:1
            
            switch simulations{sim}
                case 'DR', auv(1:vv) = AUV_Basic( sprintf('auv%d',vv), wpFile{vv});
                    
                otherwise
                    in = true(1,num_auv);
                    in(vv) = false;
                    auv(1:vv) = MR_AUV1( sprintf('auv%d',vv), bathy, wpFile{vv}, wpFile{in});
            end
        end
        
  
              
        % --- Set up Random Ploicies ----
        switch simulations{sim}
            case 'rand05', rand_policy = RandomPolicy(num_auv, sz, numcoms,  0.5); 
            case 'rand1',  rand_policy = RandomPolicy(num_auv, sz, numcoms,  1); 
            case 'rand2',  rand_policy = RandomPolicy(num_auv, sz, numcoms,  2); 
            otherwise
        end
        
        
        
        
        % --- Reset Constants ---
        speed = 3;             % ASV1 Speed                        [m/s]
        
        sim_time = 0;           % Simulation time                   [s]
        dt = 0.5;               % Simulation time step              [s]
        
        TIME(3000) = 0;         % Track Sim Time
        COMMS = false(3000,num_auv); % Track Communications
        ii = 1;
        
        turn = 1;
        finish = false;
        
        disp('Starting Simulation')
        tic
        
        % Run Simulation
        while toc < 1000
            
            % ____ Take Measurments _________________________________________
            for vv = 1:num_auv
                auv(vv) = auv(vv).Add_SensorMeasurment("Depth",  -bathy.Depth( auv(vv).State.XYZ.GT(1:2)'), sim_time);
                auv(vv) = auv(vv).Add_SensorMeasurment("Compass", auv(vv).State.XYZ.GT(3), sim_time);
            end
            
            
            % ___ Localization Update _____________________________________________
            switch simulations{sim}
                
                case 'DR' % Dear Reckoning Just Guesses
                    COMMS(ii,:) = false(1,num_auv);
                
                    
                    
                case 'PF_EKF'
                    for vv = 1: num_auv, auv(vv) = auv(vv).UpdatePF(dt); end
                    
                    COMMS(ii,:) = false(1,num_auv);
                    
                    
                    
                case 'PF_OneWay_Acoms'
                    
                    acoms = GetAcoms(auv, 1, num_auv);
                    
                    for jj = 1:num_auv
                        if jj == 1, auv(jj) = auv(jj).UpdatePF(dt);
                        else, auv(jj) = auv(jj).UpdatePF(dt, acoms{jj});
                        end
                    end
                    
                    COMMS(ii,:) = [true, false(1,num_auv-1)];
                    
                    
                    
                case 'PF_OneWay_Acoms_GT'
                    
                    auv1 = auv1.UpdatePF(dt);
                    
                    dist1 = vecnorm(auv1.State.XYZ.GT(1:2)' - auv2.State.XYZ.GT(1:2)');
                    acoms2 = {'acoms', dist1 + random('Normal', 0, 0.5), auv1.State.XYZ.GT(1:2)', [0.5, 0; 0, 0.5]};
                    auv2 = auv2.UpdatePF(dt, acoms2);
                    
                    COMMS(ii,:) = [true, false];
                    
                    
                    
                case 'PF_TwoWay_Acoms'
                    
                    acoms = GetAcoms(auv, turn, num_auv);
                    
                    for jj = 1:num_auv
                        if jj == turn, auv(jj) = auv(jj).UpdatePF(dt);
                        else, auv(jj) = auv(jj).UpdatePF(dt, acoms{jj});
                        end
                    end
                    
                    COMMS(ii,:) = false(1,num_auv);
                    COMMS(ii,turn) = true;
                    
                    turn = turn +1;
                    if turn > num_auv, turn = 1; end
                    
                    
                    
                case 'PF_ConstraintOnly_Acoms'
                    % ___ Communications Planning _________________________________________
                    for jj = num_auv: -1: 1
                        [auv(jj).Planning, turn(jj)] = auv(jj).Planning.Update_Belief(dt, auv(jj).Nav.PF);
                    end
                    
                    COMMS(ii,:) = [comms1, comms2];
                    
                    % ___ A-Comms _________________________________________________________
                    if sum(turn) ~= 1                                                        % Do not allow overlaping communications
                        for jj = 1: num_auv, auv(jj) = auv(jj).UpdatePF(dt); end
                        COMMS(ii,:) = false(1,num_auv);
                    else
                        turn = find(turn);
                        acoms = GetAcoms(auv, turn, num_auv);
                        
                        for jj = 1:num_auv
                            if jj == turn, auv(jj) = auv(jj).UpdatePF(dt);
                            else, auv(jj) = auv(jj).UpdatePF(dt, acoms{jj});
                            end
                        end
                        
                        COMMS(ii,:) = false(1,num_auv);
                        COMMS(ii,turn) = true;
                    end
                    
                 
                    
                case 'Policy'
                    turn = coms_policy(:,ii);
                    
                    if sum(turn) ~= 1                                                        % Do not allow overlaping communications
                        for jj = 1: num_auv, auv(jj) = auv(jj).UpdatePF(dt); end
                        COMMS(ii,:) = false(1,num_auv);
                    else
                        turn = find(turn);
                        acoms = GetAcoms(auv, turn, num_auv);
                        
                        for jj = 1:num_auv
                            if jj == turn, auv(jj) = auv(jj).UpdatePF(dt);
                            else, auv(jj) = auv(jj).UpdatePF(dt, acoms{jj});
                            end
                        end
                        
                        COMMS(ii,:) = false(1,num_auv);
                        COMMS(ii,turn) = true;
                    end
                    
                    
                    
                case {'rand05', 'rand1', 'rand2'}
                    turn = rand_policy(:,ii);
                    
                    if sum(turn) ~= 1                                                        % Do not allow overlaping communications
                        for jj = 1: num_auv, auv(jj) = auv(jj).UpdatePF(dt); end
                        COMMS(ii,:) = false(1,num_auv);
                    else
                        turn = find(turn);
                        acoms = GetAcoms(auv, turn, num_auv);
                        
                        for jj = 1:num_auv
                            if jj == turn, auv(jj) = auv(jj).UpdatePF(dt);
                            else, auv(jj) = auv(jj).UpdatePF(dt, acoms{jj});
                            end
                        end
                        
                        COMMS(ii,:) = false(1,num_auv);
                        COMMS(ii,turn) = true;
                    end
                    
                    
                otherwise, warning('Simulation Not Found')
            end
            
            
            
            % ___ ASV Operations __________________________________________________
            for jj = 1:num_auv
                auv(jj) = auv(jj).Navigate(speed, dt , 5);                         % Drive the AUV
                if auv(jj).Goal_Achived(0.5)||auv(jj).Nav.Done, finish = true; end  % Check for goal
            end
            
            if finish
                
                for jj = num_auv: -1: 1, stats(jj) = auv(jj).PathStates(false); end
                
                % Save Data
                for jj = 1 : ii-1
                    
%                     d{num_auv*2 + 1} = [];
                    for vv = num_auv:-1:1, d{vv} = stats(vv).err(jj); d{num_auv+1 +vv} = COMMS(jj,vv); end
                    d{num_auv+1} = TIME(jj);
                    
                    rAd = rAd.Add_Result(iter, simulations{sim}, d{:}); % Add data to the respective set
                end
                
                save(filepath, 'rAd')
                
                break 
            end
            
            
            
            % ___ Time Keeping ____________________________________________
            TIME(ii) = sim_time;
            ii = ii+1;
            sim_time = sim_time +dt;
            
            
        end
        
        
        
        fprintf('\nSimulation Time: %f [min]\nRun Time: %f [min]\n', sim_time/60, toc/60)
        disp('_____________________________________________________')
        
        rAd = rAd.ResetIndex;
        
        clear dist1 comms1 comms2 acoms1 acoms2 dt speed1 speed2 sim_time
        
        
        if plotIO
            fig = bathy.Plot_3DModel(0, 90);                                % Show AUV tracks on the Bathymetry Map
            
            for jj = 1:num_auv
                auv(jj).Nav.PlotWaypoints(fig);
                auv(jj).PlotPaths(fig);
            end
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





%%
function rand_policy = RandomPolicy(num_auv, sz, numcoms, a)

rand_policy = false(num_auv, 20000);

for vv = 1:num_auv
    rand_policy(vv, randi(sz(vv), round(numcoms(vv)*a))) = true;
end

end


function acoms = GetAcoms(auv, turn, num_auv)

for ii = num_auv: -1: 1
    dist = vecnorm(auv(ii).State.XYZ.GT(1:2)' - auv(turn).State.XYZ.GT(1:2)');
    acoms{ii} = {'acoms', dist + random('Normal', 0, 0.5), auv(turn).Nav.PF.X(1:2)', auv(turn).Nav.PF.Cov(1:2,1:2)};
end

end





