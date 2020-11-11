
% simulations = {'DR', 'PF_EKF', 'Policy', 'Full', 'eight', 'six', 'four', 'two', 'one', 'half', 'quart'};


function [rAd, coms_policy] = Simulation3(num_trials, num_auv, maxErr, simulations, bathy, saveResults, wpFile, count, plotIO)

[path, name, ext] = fileparts(saveResults);

filepath = fullfile(path, sprintf('%s_%d%s',name, count, ext));

for vv = num_auv:-1:1
    dataTypes{vv} = sprintf('Err%d',vv);
    dataTypes{num_auv+1+vv} = sprintf('Comms%d',vv);
end
dataTypes{num_auv+1} = 'Time';


rAd = RandomAssData([3000,num_trials], simulations, dataTypes);             % Instanciate RandomAssData


% ___ Communications Planning _____________________________________________
for vv = num_auv:-1:1                                                       % Create Robots for comms planning
  
   in = true(1,num_auv);
   in(vv) = false;
   auv(1:vv) = MR_AUV1( sprintf('auv%d',vv), bathy, wpFile{vv}, wpFile{in});
   
end


if plotIO
    fig = bathy.Plot_3DModel(0, 90);                                        % Show AUV tracks on the Bathymetry Map
    for jj = 1:num_auv, auv(jj).Nav.PlotWaypoints(fig); end
    drawnow
end



coms_policy = zeros(num_auv, 2000);                                         % Peralocate space for coms policy
disp('Start Communications Planning')
dt = 0.5;                                                                   % Planning time step
parfor vv = 1:num_auv                                                       % Do Planning
    
    policy = auv(vv).Planning.Plan(2000, dt, maxErr);
    policy(end+1:2000) = false;
    coms_policy(vv, :) = policy;
 
    sz(vv) = numel(policy);

end

sz = max(sz);
fprintf('\n\nFinished Communications Planning\n\n')

save(filepath, 'coms_policy', 'rAd')


% ___ Run Simulation ______________________________________________________
for iter = 1 :num_trials
    disp('____________________________________________________________________________________')
    for sim  = 1 : numel(simulations), fprintf("Iteration: %i, Simulation: %s", iter, simulations{sim})
        
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
            case 'eight', rand_policy = RandomPolicy(num_auv, sz, 8.0 /10); 
            case 'six',   rand_policy = RandomPolicy(num_auv, sz, 6.0 /10); 
            case 'four',  rand_policy = RandomPolicy(num_auv, sz, 4.0 /10);
            case 'two',   rand_policy = RandomPolicy(num_auv, sz, 2.0 /10);
            case 'one',   rand_policy = RandomPolicy(num_auv, sz, 1.0 /10);
            case 'half',  rand_policy = RandomPolicy(num_auv, sz, 0.5 /10);
            case 'quart', rand_policy = RandomPolicy(num_auv, sz, 0.25/10);
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
        
        
        % ==== Run Simulation =============================================
        tic
        while toc < 1000
            
            % ____ Take Measurments _________________________________________
            for vv = 1:num_auv
                auv(vv) = auv(vv).Add_SensorMeasurment("Depth",  -bathy.Depth( auv(vv).State.XYZ.GT(1:2)'), sim_time);
                auv(vv) = auv(vv).Add_SensorMeasurment("Compass", auv(vv).State.XYZ.GT(3), sim_time);
            end
            
            
            % ___ Localization Update _____________________________________________
            switch simulations{sim}
                
                case 'DR' % Dead Reckoning Just Guesses
                    
                    
                case 'PF_EKF'
                    for vv = 1: num_auv, auv(vv) = auv(vv).UpdatePF(dt); end
                    

                case 'Full'
                    
                    acoms = GetAcoms(auv, turn, num_auv);
                    
                    for jj = 1:num_auv
                        if jj == turn, auv(jj) = auv(jj).UpdatePF(dt);
                        else, auv(jj) = auv(jj).UpdatePF(dt, acoms{jj});
                        end
                    end
                    
                    COMMS(ii,turn) = true;
                    
                    turn = turn +1;
                    if turn > num_auv, turn = 1; end
                    
                 
                    
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
                    
                    
                    
                case {'eight', 'six', 'four', 'two', 'one', 'half', 'quart'}
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
            
            
            % ___ Check For Goal __________________________________________
            if finish
                
                for jj = num_auv: -1: 1, stats(jj) = auv(jj).PathStates(false); end
                
                % Save Data
                for jj = 1 : ii-1
                    
%                     d{num_auv*2 + 1} = [];
                    for vv = num_auv:-1:1, d{vv} = stats(vv).err(jj); d{num_auv+1 +vv} = COMMS(jj,vv); end
                    d{num_auv+1} = TIME(jj);
                    
                    rAd = rAd.Add_Result(iter, simulations{sim}, d{:}); % Add data to the respective set
                end
                
                save(filepath, 'coms_policy', 'rAd')
                
                break 
            end
            
            
            
            % ___ Time Keeping ____________________________________________
            TIME(ii) = sim_time;
            ii = ii+1;
            sim_time = sim_time +dt;
            
            
        end
        
        
        
        fprintf('  Sim Time: %.2f [min]  Run Time: %.2f [min]\n', sim_time/60, toc/60)
        
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

rAd = rAd.Eval_Stats('Time');
save(filepath, 'coms_policy', 'rAd')

fprintf('\nNumber of Communications\n')
for ii = 1: num_auv
    fprintf('AUV %d: %d\n',ii, sum(coms_policy(ii,:)))
end

% clear path
% 
% for ii = 1: num_auv
%     path(ii) = {auv(1, ii).Paths.GT};
% end
% 
% save(sprintf('/Users/jake/_OutPuts/Sim_Results/20200427_path_DR_%d.mat',num_auv), 'path')

end





%%
function rand_policy = RandomPolicy(num_auv, sz, a)

rand_policy = false(num_auv, 20000);

I = eye(num_auv);

b = round(sz * a / num_auv);

for ii = linspace(1,sz,b) 
    b = floor(ii);
    rand_policy(:,b:b+num_auv-1) = I;
end

end


function acoms = GetAcoms(auv, turn, num_auv)

for ii = num_auv: -1: 1
    dist = vecnorm(auv(ii).State.XYZ.GT(1:2)' - auv(turn).State.XYZ.GT(1:2)');
    acoms{ii} = {'acoms', dist + random('Normal', 0, 0.5), auv(turn).Nav.PF.X(1:2)', auv(turn).Nav.PF.Cov(1:2,1:2)};
end

end





