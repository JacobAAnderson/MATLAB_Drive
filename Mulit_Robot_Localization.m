 % Multi-robot localization Sim
 % Jacob Anderson
 % 1/30/2020
 
close all
clear all
clc
format long g

ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(false);                                                     % Indicate whether new data should be selescted by the user


%% Get Data Sources

num_auv = 2;

bathFile = ui.GetFile('bathymetry','mat');                                  % Get file path GUI
if isempty(bathFile), return, end                                           % End script if there isn't an input for it to use

for jj = num_auv: -1: 1
    wpFile{jj} = ui.GetFile( sprintf('waypoints%d', jj),'txt');             % Get file path GUI
end

if any(isempty(wpFile)), return, end

saveResults = ui.SaveFile('Results','mat');                                 % Get File Path to save Results


%% Create Robots
bathy = LoadBathymetry(bathFile);                                           % Load Bathymetry Map
bathy = bathy.ReduceResolution(0.5);

for vv = num_auv: -1: 1
    
%    auv(vv) = AUV_Basic( 'auv1', wpFile{vv});
   
   in = true(1,num_auv);
   in(vv) = false;
   auv(1:vv) = MR_AUV1( sprintf('auv%d',vv), bathy, wpFile{vv}, wpFile{in});
   
end

% clear bathFile wpFile vv in                                                 % Clean up workspace


%% Diplay Initial Setup
% fig = bathy.Plot_3DModel(0, 90);                                            % Show the Bathymetry Map

% for ii = 1: numel(auv) 
%     auv(ii).Nav.PlotWaypoints(fig);
%     auv(ii).PlotLocation(fig);
% end


%% Communications Planning
speed1 = 3;             % ASV1 Speed                        [m/s]

sim_time = 0;           % Simulation time                   [s]
dt = 0.5;               % Simulation time step              [s]

TIME(10000) = 0;        % Track Sim Time
COMMS = false(10000,2); % Track Communications
turn = 1;
ii = 1;

finish = false;

coms_policy = zeros(num_auv, 2000);

fprintf('\n\nStarting Communications Planning\n\n')
for jj = 1: num_auv
    fprintf('Planning for AUV on: %s\n',wpFile{jj}) 
    policy = auv(jj).Planning.Plan1(10, dt, 1.25, false);
    coms_policy(jj, 1:size(policy,2) ) = policy;
end
p = sum(coms_policy,1);

fprintf('\n\nFinished Communications Planning\n')
fprintf('Number of overlaping communications: %d\n\n', sum(p>1))


% save(sprintf('%s_CommsPolicies.mat', TODAY),'coms_policy1','coms_policy2')
% load('20200224_CommsPolicies.mat')


%% Run Simulation
disp('Starting Simulation')
tic
while toc < 1000
    
    
    % ____ Take Measurments _________________________________________
    for jj = 1:num_auv
        auv(jj) = auv(jj).Add_SensorMeasurment("Depth",  -bathy.Depth( auv(jj).State.XYZ.GT(1:2)'), sim_time);
        auv(jj) = auv(jj).Add_SensorMeasurment("Compass", auv(jj).State.XYZ.GT(3), sim_time);
    end
    

    % ___ PF Only Update __________________________________________________
%     for jj = 1: num_auv, auv(jj) = auv(jj).UpdatePF(dt); end


    % ___ Two Way Comms Update ____________________________________________
%     dist1 = vecnorm(auv(1).State.XYZ.GT(1:2)' - auv(2).State.XYZ.GT(1:2)');
%     acoms{1} = {'acoms', dist1 + random('Normal', 0, 0.5), auv(2).Nav.PF.X(1:2)', auv(2).Nav.PF.Cov(1:2,1:2)};
%     acoms{2} = {'acoms', dist1 + random('Normal', 0, 0.5), auv(1).Nav.PF.X(1:2)', auv(1).Nav.PF.Cov(1:2,1:2)};    
%     for jj = 1:num_auv, auv(jj) = auv(jj).UpdatePF(dt, acoms{jj}); end

    
    % ___ Take Turns Communicating ________________________________________
%     acoms = GetAcoms(auv, turn, num_auv);
%     
%     for jj = 1:num_auv 
%         if jj == turn, auv(jj) = auv(jj).UpdatePF(dt);
%         else, auv(jj) = auv(jj).UpdatePF(dt, acoms{jj});
%         end
%     end
%     
%     turn = turn +1;
%     if turn > num_auv, turn = 1; end


    % ___ Constraint Based Communications Planning ________________________
%      [auv1.Planning, comms1] = auv1.Planning.Update_Belief(dt, auv1.Nav.PF);
%      [auv2.Planning, comms2] = auv2.Planning.Update_Belief(dt, auv2.Nav.PF);
%      
%      COMMS(ii,:) = [comms1, comms2];


    % ___ Planning Based Comms ____________________________________________
    turn = coms_policy(:,ii);
    
    if sum(turn) ~= 1                                                        % Do not allow overlaping communications
        for jj = 1: num_auv, auv(jj) = auv(jj).UpdatePF(dt); end
    else
        turn = find(turn);
        acoms = GetAcoms(auv, turn, num_auv);
        
        for jj = 1:num_auv
            if jj == turn, auv(jj) = auv(jj).UpdatePF(dt);
            else, auv(jj) = auv(jj).UpdatePF(dt, acoms{jj});
            end
        end
    end




    % ___ ASV Operations __________________________________________________
    for jj = 1:num_auv 
        auv(jj) = auv(jj).Navigate(speed1, dt , 5);                         % Drive the AUV
        if auv(jj).Goal_Achived(0.5)||auv(jj).Nav.Done, finish = true; end  % Check for goal
    end
    
    if finish, break, end
    
    TIME(ii) = sim_time;
    ii = ii+1;
    sim_time = sim_time +dt;                                                % time keeping

    
    
    % ____ Visulization ___________________________________________________
%     if mod(sim_time, 50) == 0
%         for jj = 1: num_auv
%             auv(jj).PlotPaths(fig);
%             auv(jj).Nav.PF.PlotState(fig, 0.5);
%             auv(jj).PlotState(fig);
%             drawnow
%         end
%     end
    
%   for ii = 1: num_auv 
%       auv(jj).PlotLocation(fig);
%       drawnow
%   end
    
end

disp('_____________________________________________________')
fprintf('\nSimulation Time: %f [min]\nRun Time: %f [min]\n', sim_time/60, toc/60)
disp('_____________________________________________________')

clear dist1 comms1 comms2 acoms1 acoms2 dt speed1 speed2 sim_time


%% Plotting
fig = bathy.Plot_3DModel(0, 90);                                            % Show the Bathymetry Map

for jj = 1:num_auv
    auv(jj).Nav.PlotWaypoints(fig);
    auv(jj).PlotPaths(fig);
%   axis equal
end

%% Run Stats

figure
for jj = num_auv: -1: 1
    stats(jj) = auv(jj).PathStates(false);
    
    subplot(num_auv,1,jj)
    plot(stats(jj).err)
    title(sprintf('AUV%d',jj))
    ylabel('Error [m]')
    xlabel('time step')
end    
drawnow


%% Save Results
if isempty(saveResults), return, end                                        % End script if there isn't an input for it to use

simulations = {'DR'};
dataTypes   = {'Err1', 'Err2', 'Err3', 'Time'};

rAd = RandomAssData([ii-1,1], simulations, dataTypes);                  % Instanciate RandomAssData

for jj = 1 : ii-1
    rAd = rAd.Add_Result(1, 'DR', stats(1).err(jj), stats(2).err(jj), stats(3).err(jj), TIME(jj)); % Add data to the respective set
end
save(saveResults, 'rAd')                                                    % Save Data


figure('Name', 'Error Results', 'numbertitle', 'off')
for jj = 1: num_auv
    subplot(num_auv,1,jj)
    rAd.PlotData(simulations, sprintf('Err%d',jj), 'Time');
end




%% === FUNCTIONS ======================================================================

function acoms = GetAcoms(auv, turn, num_auv)

for ii = num_auv: -1: 1
    dist = vecnorm(auv(ii).State.XYZ.GT(1:2)' - auv(turn).State.XYZ.GT(1:2)');
    acoms{ii} = {'acoms', dist + random('Normal', 0, 0.5), auv(turn).Nav.PF.X(1:2)', auv(turn).Nav.PF.Cov(1:2,1:2)};
end

end




