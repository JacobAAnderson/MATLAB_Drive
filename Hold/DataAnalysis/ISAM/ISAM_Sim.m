%% Ecomapper ISAM with ROS
% 07/26/2019

% This Script reads in SSS data from the Ecomapper
% It Particle filters the altimiter reading on a bathymetry map
% Sends SSS images to a feature extractor node and recives resutls via ros
% Sends image patches of the identified features to a matching node via ros
% Sends odometry and landmarks to ISAM node via ros

close all
clear all
clc
format long g

% Get Data
ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(false);                                                     % Indicate whether new data should be selescted by the user

file = ui.GetFile('ISAM_Data','mat');                                       % GUI to get the name and file path of a file
if isempty(file), return, end                                               % Check if the Data Input box was cnaceled

saveFile = ui.SaveFile('ISAM_Results','mat');                               % GUI to get file path to save variable / object
if isempty(saveFile), warning('Results will Not be Saved'), end             % End script if there isn't an input for it to use


load(file)

sssData = sssData.EM_Feature_Index;                                         % Correlate reading to the Ecomapper time stamps
sssData = sssData.MakeSoanrRanges(25);                                      % make Sonar Ranges

sssData = sssData.GetSubMapIndex;

results = MakeISAMResults(bathymetry, emData);                              % Get Results Data Structure

slam = ISAM_ROS;                                                            % Start ROS


%% Run Simulation

clc

sig = 6;
count = 1;
plotting = false;
prob_threshold = 0.0;

tic
for iter  = 1: 50                                                           % {1:100}         Number of iteration for each scheem
for param = 0.9:0.9                                                         % {0.1:0.05:0.51} Sweep through Sizes of bathymetery map
for type  = 1:2                                                             % {1:3}           Compare 1) pf to 2) pf-slam to 3) slam

    fprintf("Run: %d, parm: %f, type: %d\n", iter, param, type)

    % Fresh start
    close all
    clear pf features path
    path = emData.Path.DR.utm;
%    path = emData.Path.GT.utm;

    slam.Clear_SLAMGraph;
    slam.Clear_LandMarks;
    slam.Clear_SubMap;
    slam.Clear_SLAM_Update

    slam.Publish_GPS(emData.Path.GT.utm(1,:))                               % Publish Initial GPS location

    bathymetry = bathymetry.ReduceResolution(param);                        % Down Sample Bathymetery Map

    sssData = sssData.Add_Batymetry(bathymetry);

    pf = StartPF(bathymetry, emData, 0.9, 0.2, 3);                          % Instanciate the Particle filter with 1000 particles and the Bathymetery Map

    %___Plotting___________________________________________________________
    if plotting, fig = bathymetry.Plot_3DModel(0, 90); hold on, end


% Run the Mission =======================================================================================================================
    for ii = 2: max(sssData.EM_Index(:))

        % SLAM Update Recived ---------------------------------------------
        if slam.Update.newData && type > 1
            
             %___Plotting___________________________________________________________
            if plotting, figure(fig), plot(slam.Update.pose(1), slam.Update.pose(2), '*r'); end

            %             fprintf("\nSLAM Update: Old Pose --> %f,  %f\n", pf.GetLocation)
            %             fprintf("             New Pose --> %f,  %f\n\n", slam.Update.pose(1), slam.Update.pose(2))
            %
            %             disp("SLAM COV:")
            %             disp(slam.Update.cov)
            %
            %             disp(" Covariance Determinant ")
            %             disp( det(slam.Update.cov))

            if type == 2
                pf = pf.AdjustLocation(slam.Update);
            elseif type == 3
                d_xy = slam.Update.pose - path(ii,:);
                path(ii:end,:) = path(ii:end,:) + d_xy;
            else
                disp("Too Many types")
            end

           

            slam.Clear_SLAM_Update;

        end

        if ii > max(sssData.EM_Index(:)), continue, end
        
        % Update Particle Filter ------------------------------------------
        pf = ParticleFilterUpdate(pf, emData, ii);
        if type < 3, path(ii,:) = pf.GetLocation; end

        % Look for features -----------------------------------------------
        if type > 1 && ismember(ii, [sssData.Features.EM_index])              % There is a feature at this time step

            heading = emData.FilteredData.Heading(ii);

            sssData = sssData.Feature_Make_Cov(ii, heading, emData.FilteredData.Altimeter(ii), pf.Cov(1:2,1:2));

            feature = sssData.ISAM_feature(ii, path(ii,:), pf.Cov);
            feature = sssData.Feature_CentroidfromBathy(feature, path(ii,:), heading);

            slam.AddFeature2submap(feature);

%             feature = slam.Feature_PotnetialMatch(feature, prob_threshold, false);
%             SendSLAMFeatures( slam, feature, path(ii,:), pf.Cov, fig);

 
        end


        % Align submaps -------------------------------------------------------------------------------------
        if type > 1 && ismember(ii, [sssData.SubMap.Indx])
            disp("Align Submaps")
            slam.AlignMaps(true);

            feature = slam.MergeMaps(prob_threshold, false);

            slam.Clear_SubMap;

            if ~isempty(feature(1).name)
               
                SendSLAMFeatures( slam, feature, path(ii,:), pf.Cov);
 
            end

        end

        %___Plotting___________________________________________________________
        if plotting
            figure(fig)
            p1 = plot(path(ii-1: ii,1), path(ii-1:ii,2), 'b');   % Plot Vehicle path
            p2 = plot(emData.Path.GT.utm(ii-1:ii,1), emData.Path.GT.utm(ii-1:ii,2), 'g');
%             p3 = plot(emData.Path.DR.utm(ii-1:ii,1), emData.Path.DR.utm(ii-1:ii,2), 'r');
            drawnow
        end

    end
    % =======================================================================================================================================

    %___Plotting___________________________________________________________
    if plotting, legend([p1, p2],{'SLAM','GPS'}), hold off, end

    clear speed heading dT measurement offSet d_xy in                    % clean up work space

    slam.Clear_SLAMGraph;                                                   % Clear the SLAM graph for the next run

    % Calculate path RMS error and assing to appropriate field
    switch type

        case 1
            results.RMS_err(count).PF = rmse( emData.Path.GT.utm(1:ii,:) - path(1:ii,:) );
            results.Paths(count).PF   = path(1:ii,:);
            count = count + 1;
        case 2
            results.RMS_err(count).PF_SLAM = rmse( emData.Path.GT.utm(1:ii,:) - path(1:ii,:) );
            results.Paths(count).PF_SLAM   = path(1:ii,:);

        case 3
            results.RMS_err(count).SLAM = rmse( emData.Path.GT.utm(1:ii,:) - path(1:ii,:) );
            results.Paths(count).SLAM   = path(1:ii,:);
            results.Param(count)        = param;

            count = count + 1;
    end

end

save(saveFile, 'results'  )                                                 % Save Peroidicaly incase something crashes

end
end

results = results.Truncate(count);                                          % Get Rid of extra space

save(saveFile, 'results' )                                                  % Save Final data

slam.End;

% Done !!!!
fprintf("\n\nDone!!, Simulation Time: %f\n\n\n", toc)


%% Plotting
% PlotSTATS(pf.stats, path.err, emData.timestamp );
% PlotPathes(Path_Results, Bathymetry_Map2)                                 % Plot Resultant Path
% PlotRMSerror(emData.timestamp, RMS_Error(1))



%% === Sub Routiens =======================================================================================================================


function SendSLAMFeatures( slam, feature, path, cov, fig)

if nargin > 4
    slam.Publish_Landmark( feature, path, cov, fig);
    slam.Publish_GraphOptimize(path, cov, fig);
else
    slam.Publish_Landmark( feature, path, cov);
    slam.Publish_GraphOptimize(path, cov);
end


% % Check to see if the new features have a match
% new = false(size(feature));
% for aa = 1: size(feature,2)
%     new(aa) = feature(aa).matchName == "";
% end
% 
% if any(new)  % Add feature if matching conditiona are not met
%     slam.Publish_Landmark( feature(new), path, cov, fig);
% end
% 
% if any(~new)    % Close Loop
%     slam.Publish_CloseLoop(feature(~new), path, cov, fig);
%     slam.Publish_GraphOptimize(path, cov, fig);
% end

end


% Make Sim Results data structure
function results = MakeISAMResults(bathymetry, emData)

x_range = max(bathymetry.Easting(:))  - min(bathymetry.Easting(:));
y_range = max(bathymetry.Northing(:)) - min(bathymetry.Northing(:));

results = ISAM_Results(5000);
results = results.AddMapDim(size(bathymetry.Easting), [x_range, y_range]);
results = results.AddGroundTruth( emData.Path.GT.utm);
results = results.AddDearReconing(emData.Path.DR.utm);
end


% Start Particle Filter
function pf = StartPF(bathymetry, emData, sensorNoise, vehicleNoise, compassNoise)
pf = ParticleFilter(10000, bathymetry, size(emData.RawData.timeStamp,1));  % Instanciate the Particle filter with 1000 particles ane the Bathymetery Map
pf = pf.SetLocation( [emData.Path.DR.utm(1,1), emData.Path.DR.utm(1,2)] ); % Set the initial vehilce location
pf = pf.SetNoise_Sensor('Normal',  0, sensorNoise);   % 0.9                % Set Sensor model as a probability distribuition ( "type of distribuition' , mu, sigma )
pf = pf.SetNoise_Vehicle('Normal', 0, vehicleNoise);   % 0.2               % Set vehicle transiton model as a probability distribuition ( "type of distribuition' , mu, sigma )
pf = pf.SetNoise_Compass('Normal', 0, compassNoise);                       % Set Comapss Noise
end


% Particle filter Update
function pf = ParticleFilterUpdate(pf, emData, ii)

measurement = emData.RawData.bathymetry(ii,3) + emData.RawData.tide(ii);
offSet      = emData.Path.alt.utm(ii,:) - emData.Path.GT.utm(ii,:);
speed       = emData.RawData.attitude(ii,5);
heading     = emData.FilteredData.Heading(ii);
dT          = seconds( emData.RawData.timeStamp(ii) - emData.RawData.timeStamp(ii-1) );

pf = pf.Update(speed, heading, dT, measurement, offSet);
end



%% Plotting
function Plot_Feature(fig, feature )
figure(fig)
for ii = 1: size(feature,2)
    plot(feature(ii).xyz(1), feature(ii).xyz(2),'+k')
    plotErrorEllipse(feature(ii).xyz(1:2), feature(ii).error, 0.85)
end
end


function Plot_Match(fig, feature )
figure(fig)
for ii = 1: size(feature,2)
    plot(feature.xyz(1), feature.xyz(2), '+m')
    plot(feature.matchXY(1), feature.matchXY(2), 'om')
    plotErrorEllipse(feature.xyz(1:2), feature.cov, 0.85)
end

for ii = 1: size(feature,2)
    figure
    subplot(1,2,1)
    imshow(feature(ii).image)
    
    subplot(1,2,2)
    imshow(feature(ii).matchImage)
    
end

end


% Plot Paths
function PlotPathes(paths, Bathyimetry_Map)


figure( PlotBathymetry(Bathyimetry_Map.easting, Bathyimetry_Map.northing, Bathyimetry_Map.elevation, 0, 90) );

hold on

% Plot Starting Point
plot(paths(1).GT(1,1),      paths(1).GT(1,2), 'dk');

% Plot Paths
p1 = plot(paths(1).GT(:,1),      paths(1).GT(:,2), 'g');
p2 = plot(paths(1).DR(:,1),      paths(1).DR(:,2), 'r');
p3 = plot(paths(1).PF(:,1),      paths(1).PF(:,2), 'b');
p4 = plot(paths(1).SLAM(:,1),    paths(1).SLAM(:,2), 'c');
p5 = plot(paths(1).PF_SLAM(:,1), paths(1).PF_SLAM(:,2), 'm');

% Plot End Markers
plot(paths(1).GT(end,1),      paths(1).GT(end,2), '+g');
plot(paths(1).GT(end,1),      paths(1).GT(end,2), 'og');

plot(paths(1).DR(end,1),      paths(1).DR(end,2), '+r');
plot(paths(1).DR(end,1),      paths(1).DR(end,2), 'or');

plot(paths(1).PF(end,1),      paths(1).PF(end,2), '+b');
plot(paths(1).PF(end,1),      paths(1).PF(end,2), 'ob');

plot(paths(1).SLAM(end,1),    paths(1).SLAM(end,2), '+c');
plot(paths(1).SLAM(end,1),    paths(1).SLAM(end,2), 'oc');

plot(paths(1).PF_SLAM(end,1), paths(1).PF_SLAM(end,2), '+m');
plot(paths(1).PF_SLAM(end,1), paths(1).PF_SLAM(end,2), '0m');

hold off

legend([p1,p2,p3,p4,p5],["GPS","DR", "PF", "SLAM", "PF-SLAM"])

end


% Plot Stats
function PlotSTATS(stats, err, time)

figure('Name', 'Simulation Stats', 'numbertitle','off')
subplot(4,1,1)
plot(time, err, 'r');
title( "RMS Err")

subplot(4,1,2)
plot(time, stats.mean(:,3), 'b');
title("Average Weight Values")


subplot(4,1,3)
plot(time, stats.max_min(:,1), 'b');
title("Max Weight Values")

% subplot(4,1,3)
% plot(time, stats.std(:,1), 'm');
% hold on
% plot(time, stats.std(:,2), 'c');
% hold off
% title("STD of position")

subplot(4,1,4)
plot(time, stats.std(:,3), 'g');
title("SDT of Weights")

end


% Plor RMS Error
function PlotRMSerror(time, RMS_Error)

figure('name', 'RMS Error', 'numbertitle','off');
plot(time, RMS_Error.pf, 'r')
hold on
plot(time, RMS_Error.pf_slam, 'b')
hold off
title("PF vs PF-SLAM RMS Error")
legend("PF", "PF-SLAM")

end


