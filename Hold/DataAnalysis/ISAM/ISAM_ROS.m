% ISAM_ROS data class
% Jake Anderson
% 10/3/2019

classdef ISAM_ROS < handle
    
    properties
        Pub
        Sub
        Update
        LandMarks
        SubMap
    end
    
    methods
        
        % Constructor
        function obj = ISAM_ROS(name)
          
            rosinit('Me');                                                  % Conect to ros master
            
            % Initialize Publishers and subscriber
            if nargin < 1
                obj.Pub.odom   = rospublisher('/odom',          'nav_msgs/Odometry');            % Publish odometry Data
                obj.Pub.gps_lm = rospublisher('/gps_landmark',  'nav_msgs/Odometry');            % Publish SSS raw image
                obj.Pub.opt_gh = rospublisher('/optimizeGraph', 'nav_msgs/Odometry');            % Publish odometry Data
                obj.Pub.sss_lm = rospublisher('/sss_landmark',  'sensor_msgs/PointCloud');       % Publish SSS identified features
                obj.Pub.matchs = rospublisher('/sss_match',     'sensor_msgs/PointCloud');       % Publish SSS Feature matches
                obj.Pub.clg    = rospublisher('/clearGraph',    'std_msgs/String');              % Clear SLAM Graph
                
                obj.Sub.slam_poses = rossubscriber('/slam_update_pose', @obj.SLAM_Pose_Update); % Subscribe to a SSS image feature identification mesages
                
            else
                obj.Pub.odom   = rospublisher([name,'/odom'],          'nav_msgs/Odometry');            % Publish odometry Data
                obj.Pub.gps_lm = rospublisher([name,'/gps_landmark'],  'nav_msgs/Odometry');            % Publish SSS raw image
                obj.Pub.opt_gh = rospublisher([name,'/optimizeGraph'], 'nav_msgs/Odometry');            % Publish odometry Data
                obj.Pub.sss_lm = rospublisher([name,'/sss_landmark'],  'sensor_msgs/PointCloud');       % Publish SSS identified features
                obj.Pub.matchs = rospublisher([name,'/sss_match'],     'sensor_msgs/PointCloud');       % Publish SSS Feature matches
                obj.Pub.clg    = rospublisher([name,'/clearGraph'],    'std_msgs/String');              % Clear SLAM Graph
                
                obj.Sub.slam_poses = rossubscriber('/slam_update_pose', @obj.SLAM_Pose_Update); % Subscribe to a SSS image feature identification mesages
                
            end
            
            % Initialize Data structures
            obj.Update.newData = false;
            obj.Update.pose = [ 0, 0 ];
            obj.Update.cov = zeros( 3, 3);
            obj.LandMarks = sssFeature;
            pause(2)                                                        % Wait to make sure publisher is registered
            
        end
        
        
        % Send GPS Location
        function Publish_GPS(obj, gps)
            
            odomMsg = rosmessage(obj.Pub.gps_lm);
            
            odomMsg.Pose.Pose.Position.X = gps(1);
            odomMsg.Pose.Pose.Position.Y = gps(2);
            odomMsg.Pose.Pose.Position.Z = 0;
            odomMsg.Pose.Covariance(1) = 0;
            odomMsg.Pose.Covariance(2) = 0;
            
            send(obj.Pub.gps_lm, odomMsg)
            
        end
        
        
        % Send Feature message
        function Publish_Landmark(obj, features, vehilcePose, err, fig)
            
            if isempty(obj.LandMarks(1).name)
                c = 1;
            else
                c = max(size(obj.LandMarks)) + 1;
            end
            
            ids = c: c + size(features,2) - 1;
            
            obj.LandMarks(ids) = features;
            
            fprintf("Send Feature : %d\n", ids)
            
            odomMsg    = rosmessage(obj.Pub.odom);
            sssMsg     = rosmessage(obj.Pub.sss_lm);
            
            for ii = 1: size(features,2)
                
                pointMsg   = rosmessage('geometry_msgs/Point32');
                channelMsg = rosmessage('sensor_msgs/ChannelFloat32');
                
                pointMsg.X = features(ii).xyz(1);
                pointMsg.Y = features(ii).xyz(2);
                pointMsg.Z = features(ii).xyz(3);
                
                channelMsg.Values(1) = ids(ii);
                channelMsg.Values(2:5) = features(ii).cov(:);
                
                sssMsg.Channels(ii) = channelMsg;
                sssMsg.Points(ii) = pointMsg;
                
                odomMsg.Pose.Pose.Position.X = features(ii).vehicle(1);
                odomMsg.Pose.Pose.Position.Y = features(ii).vehicle(2);
                odomMsg.Pose.Pose.Position.Z = 0;
                odomMsg.Pose.Covariance(1:9) = features(ii).vehicle_cov(:);
                
                send(obj.Pub.odom, odomMsg)
                send(obj.Pub.sss_lm, sssMsg)
                
                if nargin > 4
                    figure(fig)
                    plot(features(ii).xyz(1), features(ii).xyz(2),'+k')
                    plotErrorEllipse(features(ii).xyz(1:2), features(ii).error, 0.85,'k')
                    drawnow
                end
                
                pause(0.1)
                
            end
        end
        
        
        % Send Loop closures
        function Publish_CloseLoop(obj, features, vehilcePose, err, fig)
            
            odomMsg  = rosmessage(obj.Pub.odom);
            matchMsg = rosmessage(obj.Pub.matchs);
            
            for ii = 1: max(size(features))
                
                ids = find( ismember([obj.LandMarks.name], features(ii).matchName) );
                
                obj.LandMarks(ids).ind = features(ii).ind;                  % Update When this Feature was last seen
                
                if isempty(ids)
                    warning("Not enought Features")
                    return
                end
                
                fprintf("Close Loop with Feature: %d\n", ids)
                
                
                pointMsg   = rosmessage('geometry_msgs/Point32');
                channelMsg = rosmessage('sensor_msgs/ChannelFloat32');
                
                pointMsg.X = features(ii).xyz(1);
                pointMsg.Y = features(ii).xyz(2);
                pointMsg.Z = features(ii).xyz(3);
                
                channelMsg.Values(1) = ids;
                channelMsg.Values(2:5) = features(ii).cov(:);
                
                matchMsg.Channels(ii) = channelMsg;
                matchMsg.Points(ii) = pointMsg;
                
                
                
                
                odomMsg.Pose.Pose.Position.X = features(ii).vehicle(1);
                odomMsg.Pose.Pose.Position.Y = features(ii).vehicle(2);
                odomMsg.Pose.Pose.Position.Z = 0;
                odomMsg.Pose.Covariance(1:9) = features(ii).vehicle_cov(:);
                
                send(obj.Pub.odom, odomMsg)
                send(obj.Pub.matchs, matchMsg)
                
                if nargin > 4
                    figure(fig)
                    plot(features(ii).xyz(1), features(ii).xyz(2),'+m')
                    plotErrorEllipse(features(ii).xyz(1:2), features(ii).error, 0.85)
                    
                    plot(obj.LandMarks(ids).xyz(1), obj.LandMarks(ids).xyz(2), '+r')
                    plotErrorEllipse(obj.LandMarks(ids).xyz(1:2), obj.LandMarks(ids).error, 0.85)
                end
                
                
                pause(1)    % Allow time for batch optimization
                
            end
        end
        
        
        % Send Feature message
        function Publish_GraphOptimize(obj, vehilcePose, err, fig)
            
            disp("Optimize Graph")
            
            odomMsg = rosmessage(obj.Pub.opt_gh);
            
            odomMsg.Pose.Pose.Position.X = vehilcePose(1);
            odomMsg.Pose.Pose.Position.Y = vehilcePose(2);
            odomMsg.Pose.Pose.Position.Z = 0;
            odomMsg.Pose.Covariance(1:9) = err(:);
            
            send(obj.Pub.opt_gh, odomMsg)
            
            if nargin > 3
                figure(fig)
                plot(vehilcePose(1), vehilcePose(2),'dw')
                plotErrorEllipse(vehilcePose, err(1:2,1:2), 0.85, 'w')
            end
            
            pause(1)
            
        end
        
        
        % SLAM Pose Update
        function SLAM_Pose_Update(obj, ~, msgs)
            
            ids = [obj.LandMarks.ind];
            
            obj.Update.newData = true;
            obj.Update.pose = [ msgs.Pose.Pose.Position.X, msgs.Pose.Pose.Position.Y ];
            obj.Update.cov = reshape( msgs.Pose.Covariance(1:9), 3, 3);
            obj.Update.ind = max(ids);
%             assignin('base', 'slamUpdate', slamUpdate);
            
        end
        
        
        % SLAM Feature Update ??
        function SLAM_Feature_Update(susbscriber, msgs)
            
            disp(susbscriber(:))
            disp(msgs(:))
            
        end
        
        
        % Feature matches
        function feature = Feature_PotnetialMatch(obj, feature, prob_threshold, plot_io)
            
            % Eliminate potencial matches that have havent been seen yet or were just seen
            for jj = 1: max(size(feature))
                
                in = [obj.LandMarks.times] < feature(jj).times & abs([obj.LandMarks.times] - feature(jj).times) > seconds(20);
                
                if ~any(in), continue, end
                
                proposed = obj.LandMarks(in);
                
                
                if size(proposed,2) > 9
%                     ass = ICP_Matching(feature,  proposed);
                end
                
                
                indx = Covariance_mathcing(feature(jj), proposed, plot_io);
                
                if ~any(indx), continue, end
                
                proposed = proposed(indx);
                
                
                
%                 probs = obj.GetMatchingProb(feature(jj).name, [proposed.ImageName])
%                 in = probs >= prob_threshold;

                for ii = 1 : max(size(proposed))
                    probs(ii) = prob_from_cov(feature(jj).xyz(1:2), feature(jj).error, proposed(ii).xyz(1:2), proposed(ii).error); 
                end
                in = probs == max(probs) & probs > prob_threshold;
                
                
                if ~any(in), continue, end
                if sum(in) > 2, pause, end
                
                try
                feature(jj).matchXY    = proposed(in).xyz;
                feature(jj).matchName  = proposed(in).name;
                feature(jj).matchTime  = proposed(in).times;
                feature(jj).matchImage = proposed(in).image;
                catch
                    warning(" ISAM Matches??")
                    pause
                end
                
            end
        end
        
        
        % Add Feature to submap
        function AddFeature2submap(obj, feature)
            
            if isempty(obj.SubMap(1).name)
                c = 1;
            else
                c = max(size((obj.SubMap))) + 1;
            end
            
            for ii = 1: max(size(feature))
                obj.SubMap(c+ii-1) = feature(ii);
            end
            
        end
                
        
        % Mereg Submap into the Global Map
        function AlignMaps(obj, plot_io)
            
            % Check for empty maps
            if isempty(obj.LandMarks(1).name) || isempty(obj.SubMap(1).xyz), return, end
            
            new_xyz    = reshape([obj.SubMap.xyz],    3, []);
            global_xyz = reshape([obj.LandMarks.xyz], 3, []);
            
            % Check that there are enough points for transform
            if size(new_xyz,2) < 3 || size(global_xyz,2) < 3, return, end
            
            [TT, err] = TR_ICP(global_xyz, new_xyz);                        % Do translation only ICP
            
            if plot_io
                fprintf('\nER: %f\n',min(err))
                disp('TT')
                disp(TT(err == min(err)))
            end
            
            % Check that alignments is good
            if min(err) > 6, return, end
            
            if plot_io, disp("Small Errors"), end
            
            tt = [TT(err == min(err),:), 0];                                % Apply transform to th enew points
            icp_xyz = new_xyz - repmat(tt', 1, size(new_xyz,2));
            
            for ii = 1: size(new_xyz,2)
                obj.SubMap(ii).xyz = icp_xyz(:,ii)';
            end
            
            if plot_io
                figure
                plot3(global_xyz(1,:), global_xyz(2,:), global_xyz(3,:),'og')
                hold on
                plot3(new_xyz(1,:), new_xyz(2,:), new_xyz(3,:),'*r')
                plot3(icp_xyz(1,:), icp_xyz(2,:), icp_xyz(3,:),'*b')
                hold off
                view(0,90)
                drawnow
            end
        end
        
        
        % Merge Submaps and return features
         function features = MergeMaps(obj, prob_threshold, plot_io)
            
            features = obj.SubMap;
            
            if isempty(obj.SubMap(1).name) || isempty(obj.LandMarks(1).name), return, end
            
            % Eliminate potencial matches that have havent been seen yet or were just seen
            for jj = 1: max(size(features))
                                
                indx = Covariance_mathcing(features(jj), obj.LandMarks, plot_io);
                
                if ~any(indx), continue, end
                
                proposed = obj.LandMarks(indx);
              
                probs = [];
                for ii = 1 : size(proposed,2)
                    probs(ii) = prob_from_cov(features(jj).xyz(1:2), features(jj).error, proposed(ii).xyz(1:2), proposed(ii).error); 
                end
                in = probs == max(probs) & probs > prob_threshold;
                
                
                if ~any(in), continue, end
                if sum(in) > 2, pause, end
                
                features(jj).matchXY    = proposed(in).xyz;
                features(jj).matchName  = proposed(in).name;
                features(jj).matchTime  = proposed(in).times;
                features(jj).matchImage = proposed(in).image;

            end
            
            
        end
        
        
        % Reset the SLAM data ---------------------------------------------
        % Send Clear SLAM Graph Message
        function Clear_SLAMGraph(obj)
            
            msg = rosmessage(obj.Pub.clg);
            msg.Data = "Dam the torpidos!!";
            send(obj.Pub.clg, msg)
            pause(1)
            
            obj.LandMarks = [];
            
        end
        
        
        % Clear SLAM Features
        function Clear_LandMarks(obj)
            obj.LandMarks = [];
            obj.LandMarks = sssFeature;
        end
        
        
        % Clear SLAM Update
        function Clear_SLAM_Update(obj)
            
            obj.Update.newData = false;
            obj.Update.pose = [ 0, 0 ];
            obj.Update.cov = zeros( 3, 3);
            
        end
        
        
        % Remove all features from submap
        function Clear_SubMap(obj)
            
            obj.SubMap = [];
            
            obj.SubMap = sssFeature;
            
        end
        
        
    end
    
    methods (Static)
        function End
            rosshutdown                                                     % Shutdown ROS
        end
        
    end
    
    
end


% Gepmetric matching Probabilty
function probs = prob_from_cov(mu1, cov1, mu2, cov2)

probs = (mvnpdf(mu1,mu2,cov2)/mvnpdf(mu2,mu2,cov2)) * (mvnpdf(mu2,mu1,cov1)/mvnpdf(mu1,mu2,cov2));

if probs > 0.00001, disp(probs), end

end


% Find overlapping covariance ellipses
function indx = Covariance_mathcing(feature, proposed, plot_io)

indx = false(1, max(size(proposed)));

a1 = ErrorEllipse(feature.xyz(1:2), 1.3 * feature.error, 5);

if plot_io
    fig = figure;
    plot(a1(1,:), a1(2,:), 'k')
    hold on
end

for ii = 1 : max(size(proposed))
    
    a2 = ErrorEllipse(proposed(ii).xyz(1:2), 1.3 * proposed(ii).error, 5);
    in = inpolygon( a1(1,:), a1(2,:), a2(1,:), a2(2,:) );
    
    if any(in), indx(ii) = true; end
    
    if plot_io  && any(in)
        plot(a2(1,:), a2(2,:), 'g')
    elseif plot_io  && ~any(in)
        plot(a2(1,:), a2(2,:), 'r')
    end
    
end

if plot_io
    pause(0.5)
    close(fig)
end

end


% Evaluate Proposed mathces with Iterative Closest Point
function ind = ICP_Matching(new, old)
ind = [];

new_xyz = [reshape([new.xyz],3,[]), reshape([old(end-3:end).xyz],3,[])];
test_xyz = reshape([old(1:end-4).xyz],3,[]);

new_xyz(3,:) = 0;
test_xyz(3,:) = 0;

[TR, TT, ER] = icp(new_xyz, test_xyz);
                
fprintf("\nMin ICP Err: %f\n TR:\n", min(ER) )
disp(TR)
disp("TT:")
disp(TT)

end

