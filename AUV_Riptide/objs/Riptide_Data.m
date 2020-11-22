% Riptide Data
% Jacob Anderson
% Robotic Decition Making Laboratory (RDML)
% August 20, 2020

classdef Riptide_Data
    
    properties
        acomms
        course
        commsPlanner
        data
        filteredData
        manifest
        name
        orical
        path
        pAcommsHandler        
        rawAcomms
        rawCourse
        rawData
        tbn
        vehicleModel
        states
    end
    
    methods
        
        % Instanciate Class by loading Raw Data
        function obj = Riptide_Data(file, offset_file, name_, varargin), warning on backtrace
            
            fun = @(file_) Riptide_DataLog2Mat(file_, offset_file, varargin{:});
            obj.rawData = GetData(file, fun, 'csv');
            obj.data = obj.rawData;
            
            if nargin >= 2, obj.name = name_; 
            else, obj.name = 'Riptide';
            end
            
        end
        
        
        % Read in Acoustic Communications data
        function obj = Add_Acomms(obj,file, varargin)
            
            fun = @(file_) Riptide_AcommsLog2Mat(file_, varargin{:});
            obj.rawAcomms = GetData(file, fun, 'csv');
            obj.acomms = obj.rawAcomms;
            
        end
        
        
        function obj = Add_CommsPlanner(obj)
            
            obj.commsPlanner = CommunicationPlanning(obj.name);
            
        end
        
        
        function obj = Add_CommsPloicy(obj, name,  policy, dt)
            
            acommsTimes = obj.acomms.recived_msg;
            
            obj.acomms.policy.(name) = false( size(acommsTimes) );
            
            % -- Policy is a datetime array --
            if isa( policy, 'datetime')
                
                threshole = seconds(dt);
                
                for time = policy
                    
                    dt = abs(time - acommsTimes);
                    
                    idx = dt <= threshole & dt == min(dt);
                   
                    obj.acomms.policy.(name)(idx) = true;
                   
               end
                
                
            % -- Policy is a boolean array -- 
            elseif isa( policy, 'logical')
                
                if numel(policy) == 1
                    obj.acomms.policy.(name)(:) = policy;
                    
                else 
                    obj.acomms.policy.(name) = policy;
                    
                end
            
            % -- Policy is the wrong data type --
            else
                warning("Communication Policy is of an un-accepted type")
                disp('Communications palicy is false')
                
            end
           
        end
        
        
        function obj = Add_OWTOF(obj, owtof)
            dtfs = DateTime_Funs;
            
            for tof = owtof
                
                msgs_tof = tof.msgs;
                
                for ii = 1: numel(obj.rawAcomms)
                    
                    for member = {'recived_msg', 'sent_msg'}
                        
                        msgs_rA = obj.rawAcomms(ii).(member{1});
                        
                        [Lia,Locb] = dtfs.Ismember(msgs_tof, msgs_rA);
                        
                        if ~any(Lia),continue, end
                        
                        Locb(Locb == 0) = [];
                        
                        dtfs.SetError(msgs_tof(Lia), msgs_rA(Locb))
                        
                        obj.rawAcomms(ii).tof(Locb,:) = tof.time(Lia);
                        obj.rawAcomms(ii).gt(Locb,:) = tof.gt(Lia);
                    end
                    
                end
            end
            
        end
        
        
        % Read in pAcommsHandler debug log
        function obj = Add_pAcommsHandler(obj, file)
            
            dtfs = DateTime_Funs;
            
            % --- Read in logs ---
            fun = @(file_) Riptide_pAcommsHandler2Mat(file_);
            obj.pAcommsHandler = GetData(file, fun, 'txt');
 
            
            % --- Match data to Acoms logs ---
            if isempty(obj.rawAcomms), return, end                          % End function if acomms logs are not present
            
            fprintf("\nMatching pAcomms Logs\n")
            
            for log = obj.pAcommsHandler'                                   % Iterate through pAcomm data
                
                fprintf('\nProcessing log: %s\n', log.log)
                
                if isempty(log.msg) 
                    fprintf("\tlog is empty\n")
                    continue 
                end
                
                msgs_pAcom = log.msg.msg;                                   % pAcomm messages
                
                fprintf("\tNumber of pAcommsHandler messages: %d\n", numel(msgs_pAcom))
                
                for ii = 1: numel(obj.rawAcomms)                            % Iterate through Acomms log data
                    
                    for member = {'recived_msg', 'sent_msg';
                                  'rx_time',     'tx_time' }
                              
                        msgs_rA = obj.rawAcomms(ii).(member{1});            % Acomms log data
                        
                        [Lia,Locb] = dtfs.Ismember(msgs_pAcom, msgs_rA);    % See if any of the data matches
                        
                        if ~any(Lia), continue, end
                        
                        Locb(Locb == 0) = [];
                        
                        dtfs.SetError( msgs_pAcom(Lia), msgs_rA(Locb) )
                        
                        obj.rawAcomms(ii).pAcomms_tof(Locb,:) = log.msg.(member{2})(Lia,:);
                        
                        
                        msgs_rA(isnat(msgs_rA)) = [];
                        
                        fprintf("\tNumber of Acomms %s messages: %d\n", member{1}, numel(msgs_rA))
                        fprintf("\t%d of %d messages mached in %s\n", sum(Lia), numel(msgs_rA), obj.rawAcomms(ii).log)
                        
                    end
                    
                end
                
            end
            
        end
        
        
        % Add TBN and other patths
        function obj = Add_Path(obj, path_, type, coord_syst)
            
            path_( any(isnan(path_), 2), :) = [];
            path_( any(path_ == 0,   2), :) = [];
            
            
            switch lower(coord_syst)
                
                case {'lat-lon', 'latlon'}
                    
                    lat = path_(:,1);
                    lon = path_(:,2);
                    
                    [xx, yy] = deg2utm(lat , lon);
                    
                    
                case 'utm'
                    
                    xx = path_(:,1);
                    yy = path_(:,2);
                    
                    utmzone = obj.path.GT.utmzone;
                    utmzone = repmat(utmzone(1,:), size(xx) );
                    
                    [lat, lon] = utm2deg(xx, yy, utmzone);
                    
                    
                otherwise
                    warning('Un-reckognized coordinate system')
                    
            end
            
            obj.path.(type).lat_lon = [lat, lon];
            obj.path.(type).utm     = [xx,yy];
            
        end
        
        
        % Read in mission files to get the waypoints
        function obj = Add_Waypoints(obj, file)
            
            fun = @(file_) GetWaypoints(file_);
            obj.rawCourse = GetData(file, fun, 'txt');
            obj.course = obj.rawCourse;
            
        end
        
        
        % Provide the vehicle with a particle filter for Terrain Based Navigation annalisys
        function obj = Add_TBN(obj, map, didel)
            
            if nargin == 2
                didel = 1;
            end
            
            obj.tbn = DEC_TBN(obj.name, 1000 );                              % Instanciate Paricle Filter
            
            obj.tbn = obj.tbn.Add_Map(map);                                   % Provide Bathymetry map
            obj.tbn = obj.tbn.SetLocation(obj.path.GT.utm(1,:), [3,0;0,3]);   % Tell it where the vehicle is starting from
            

            % --- Set vehicle transiton model as a probability distribuition ( "type of distribuition' , mu, sigma ) ---
            for prop = {'acoms',  'altimeter', 'compass' 'speed';           % --> Type of noise
                        'Normal', 'Normal',    'Normal', 'Normal';          % --> Default distribion type
                         0,        0,           0,        0;                % --> Default mu
                         0.5,      2.5,         30,       2.0}              % --> Default sigma
                     
                     
%                      if  isfield( obj.vehicleModel, prop{1} )
%                          mu     = didel * obj.vehicleModel.(prop{1}).mu;
%                          sig    = didel * obj.vehicleModel.(prop{1}).sigma;
%                          name_  = obj.vehicleModel.(prop{1}).DistributionName;
%                          obj.tbn = obj.tbn.SetNoise(prop{1}, name_, mu, sig);
%                          
%                          fprintf("%s: mu: %f,  sig: %f\n", prop{1}, mu, sig)
%                          
%                      else
%                          obj.tbn = obj.tbn.SetNoise(prop{1}, prop{2}, prop{3}, prop{4});
%                      end
                     
                     obj.tbn = obj.tbn.SetNoise(prop{1}, prop{2}, prop{3}, prop{4});
                     
            end
           
            % Create perfect altimiter profile for testing
            
            gt = obj.path.GT.utm;
            dr = obj.path.DR_Corrected.utm;
            
            obj.orical.gtAlt = map.Depth(gt);
            obj.orical.drAlt = map.Depth(dr);
            
            
            
        end
        
        
        function obj = CorrectCourse(obj, geotiff)
            
            if ~isfield(obj.vehicleModel, 'compassFit_inv') || ~isfield(obj.vehicleModel, 'speedFit')
                warning("Compass model and/or speed are not avalalbe")
                return
            end
            
            path_  = obj.course.waypoints;                                  % Path needs to be lon-lat
            speed = obj.vehicleModel.speedFit;
            
            [xx, yy, zone] = deg2utm(path_(1,2), path_(1,1));               % Get Starting location in utms
            
            [distance, heading, ~] = WayPoint2DeadReckoning(path_, speed, 'm/s');
            
            f = obj.vehicleModel.compassFit_inv;
            
            heading = feval(f,heading);
            
            x = distance' .* [sind(heading), cosd(heading)];
            
            x = cumsum([xx,yy; x]);
            
            [lat, lon] = utm2deg( x(:,1), x(:,2), repmat(zone, size(x,1), 1) );
            
            obj.course.waypoints = [lon, lat];
            
            
            % Show transform
            if nargin == 2
                
                gt = obj.path.GT.lat_lon;
                
                figure('name', sprintf('%s Course Correction', obj.name) )
                geotiff.Show;
                hold on
                p1 = plot(gt(:,2), gt(:,1));
                p2 = plot(path_(:,1), path_(:,2));
                p3 = plot(lon, lat);
                hold off
                legend([p1,p2,p3], ["GT", "Origonal Course", "Corrected Course"])
                drawnow
                
            end
            
            
        end
        
        
        % Show Mission Manifest in the command Window
        function Disp_Manifest(obj)
            
            bar = "___________________________________________________________________________________________________________";
            fprintf('\n\n%s\n',bar)
            fprintf('%s Mission Manifest\n', obj.name)
            disp(obj.manifest)
            fprintf('%s\n\n',bar)
            
        end
      

        % Kalman filter altimiter data
        function [alt, obj] = Filter_Altimiter(obj, ii)
            
            % --- Filter position ----
            speed   = obj.data(1).attitude(ii,5);
            heading = obj.data(1).attitude(ii,3);
            dt      = obj.data.timeStep(ii);
            
            
            z = obj.tbn.X(1:2);                                  % Estimated location of the vehicle
            R = obj.tbn.cov(1:2,1:2);
            
            Q=eye(2);                                       % covariance of process
            
            f=@(x)[x(1) + speed * dt * sind( heading );      % nonlinear state equations
                   x(2) + speed * dt * cosd( heading )];  
            
            h=@(x) x;                                       % measurement equation
            
            x = obj.states.x;                               % state
            P = obj.states.P;                               % state covariance 
            
            
            [x, P] = ekf(f,x,P,h,z,Q,R);                    % ekf
            
            obj.states.x = x;
            obj.states.P = P;
            
%             obj.path.EKF.utm(ii,:) = x';
            
            
            
            % ---- Filter Altitude ----
            
            x = obj.tbn.Map.Depth(x');          % Altitude Estimat
            
            z = obj.filteredData.altitude(ii);      % Altimiter Measurement
            
            alt = (x + z)/2;
            
            obj.filteredData.altitude(ii) = alt;
            
            return
            
            
            sig = obj.filteredData.altSig;          % Error in the state
            
            A = 1;                                  % Process Model
            B = 1;                                  % Control Model
            C = 1;                                  % Observation Matrix
            
            Q = 1;                                  % Uncertainty in the process model
            R = obj.vehicleModel.altimeter.sigma;   % Measurment Standard Deviation
            
            u = 0;                                  % Control input
                    
            
            [alt, sig] = Kalman_Filter(x, sig, u, z, A, B, C, Q, R);        % Filter the measurment
            
            obj.filteredData.altSig = sig;                                  % Keep track of the error
            
            obj.filteredData.altitude(ii) = alt;
            
        end
        
        
        % Get Vehicle Path in UTM an degrees lat-lon
        function obj = Get_VehiclePaths(obj)
            
            % -- Clear old paths --
            obj.path = [];
            
            
            % ----- Ground true paths ------
            gps0 = obj.data(1).vehicle(:,1:2);                              % Get Vehicle Location
           
            
            % Smooth out the GPS data to determin the true vehicle heading
            gps1 = zeros(size(gps0));
            n = size(gps0,1);
            sw = SlidingWindow(20, n);
            
            for ii = 1:n
                
                data_ = gps0(sw.D1(ii),:);
                gps1(ii,:) = mean(data_);
                
            end
            
            [xx, yy, utmzone] = deg2utm( gps1(:,1), gps1(:,2) );
            
            obj.path.GT.lat_lon = gps1;
            obj.path.GT.utm     = [xx,yy];
            obj.path.GT.utmzone = utmzone;
            
            
            % ----- Dead Rekoning Paths -----
            speed     = obj.data(1).attitude(:,5);
            heading   = obj.data(1).attitude(:,3);
            timestamp = obj.data(1).timeStamp;
            
            timeStep_ = [0; seconds(timestamp(2:end) - timestamp(1: end-1))];
            
            deltaX_dr = speed .* timeStep_ .* sind( heading );
            deltaY_dr = speed .* timeStep_ .* cosd( heading );
            
            x = cumsum(deltaX_dr);
            y = cumsum(deltaY_dr);
            
            
            
            [dr_lat, dr_lon] = utm2deg(x + xx(1), y+yy(1), utmzone);
            
            obj.path.DR.utm(:,1) = x + xx(1);
            obj.path.DR.utm(:,2) = y + yy(1);
            obj.path.DR.lat_lon  = [dr_lat, dr_lon];
            
            obj.data.timeStep = timeStep_;
            
            obj.states.x = obj.path.GT.utm(1,:);
            obj.states.P = eye(2) * 3;
            
            % -- tracke EKF path --
            obj.path.EKF.utm     = NaN( size( obj.path.GT.utm) );
            obj.path.EKF.lat_lon = NaN( size( obj.path.GT.utm) );
            
            obj.path.EKF.utm(1,:)     = obj.path.GT.utm(1,:);
            obj.path.EKF.lat_lon(1,:) = obj.path.GT.utm(1,:);
            
            
            % -- Do corrected DR ?? --
            if ~ ( isfield( obj.filteredData, 'speed') || isfield( obj.filteredData, 'heading') ), return, end
            
            if isfield( obj.vehicleModel, 'speed'),   speed   = obj.filteredData.speed;   end
            if isfield( obj.filteredData, 'heading'), heading = obj.filteredData.heading; end
            
            deltaX_dr = speed .* timeStep_ .* sind( heading );
            deltaY_dr = speed .* timeStep_ .* cosd( heading );
            
            x = cumsum(deltaX_dr);
            y = cumsum(deltaY_dr);
            
            [xx, yy, utmzone] = deg2utm( obj.data(1).vehicle(:,1), obj.data(1).vehicle(:,2) );
            
            [dr_lat, dr_lon] = utm2deg(x + xx(1), y+yy(1), utmzone);
            
            obj.path.DR_Corrected.utm(:,1) = x + xx(1);
            obj.path.DR_Corrected.utm(:,2) = y + yy(1);
            obj.path.DR_Corrected.lat_lon  = [dr_lat, dr_lon];
            
            
            
        end
        
        
        % Retrive the path recorded by the particle filter
        function obj = Get_TBNpath(obj, type)
            
            utm = obj.tbn.path(1:obj.tbn.path_idx-1,1:2);
            
            zone = obj.tbn.Map.UTM_Zone;
            
            [lat, lon]  = utm2deg( utm(:,1), utm(:,2), repmat(zone, size(utm,1),1 )) ;
            
            if nargin < 2
                obj.path.tbn.utm = utm;
                obj.path.tbn.lat_lon = [lat, lon];
                
            else
                obj.path.(type).utm = utm;
                obj.path.(type).lat_lon = [lat, lon];
                
            end
            
        end
        
        
        % Get DR Path Structure
        function obj = Get_AltimeterPath(obj)
            
            for ii = 1: size(obj.rawData,2)
                
                [xx, yy, utmzone] = deg2utm(obj.rawData.bathymetry(:,1), obj.rawData.bathymetry(:,2));
                
                path_diff = obj.path(ii).DR.utm - obj.path(ii).GT.utm;
                
                obj.path(ii).alt.lat_lon = obj.rawData.bathymetry(:,1:2);
                obj.path(ii).alt.utm     = [xx, yy];
                obj.path(ii).alt.utmzone = utmzone;
                obj.path(ii).alt.DR      = obj.path(ii).alt.utm + path_diff;
                
            end
        end
        
        
        % Create a mission manifest of the data that has been read in
        function obj = Get_Manifest(obj, varargin)
 
            dispIO    = false;
            skipIdel  = false;
            alt_acoms = false;
            
            for var = varargin
                switch lower(var{1})
                    case "dispio",       dispIO    = true;
                    case "skipidel",     skipIdel  = true;
                    case "alt & acomms", alt_acoms = true;
                    case "all"
                        dispIO    = true;
                        skipIdel  = true;
                        alt_acoms = true;
                    otherwise
                        fprintf('\n\nRiptide Data -> Get_Manifest: Unknown Input " %s "\n\n', var{1})
                end
            end
            
            
            T_d = missionTabe(obj.rawData,   skipIdel);
            T_a = missionTabe(obj.rawAcomms, skipIdel);
            
            in = false(size(T_a,1),1);
            
            sz = [100,7];
            varTypes = {'string', 'datetime',  'datetime', 'string', 'datetime',  'logical', 'single'};
            varNames = {'Mission','Start Time','End Time', 'Type',   'Mid-Point', 'Has TOF', '#'};
            
            T = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
            
            % Determin which data logs have an accompaning Acomms Log
            count = 1;
            for indx = 1:size(T_d,1)
                
                miss = T_d.('Mission')(indx);
                time = T_d.('Mid-Point')(indx);
                
                type = "Alt";
                hasTOF = false;
                for jj = 1:size(T_a,1)
                   
                    if strcmp(miss, T_a.Mission(jj)) && T_a.('Start Time')(jj) < time && time < T_a.('End Time')(jj)
                        type = "Alt & Acomms";
                        hasTOF = T_a.('Has TOF')(jj);
                        in(jj) = true;
                        break
                    end
                end
                
                
                if ~alt_acoms || (alt_acoms && type == "Alt & Acomms")
                    T(count,:) = {miss, T_d.('Start Time')(indx), T_d.('End Time')(indx), type, time, hasTOF, count};
                    count = count + 1;
                end
            end
            
            % Collect the Acomms Logs that are not associated with a data log
            if ~alt_acoms
                for in = find(~in)'
                    T(count,:) = {T_a.Mission(in), T_a.('Start Time')(in), T_a.('End Time')(in), "Acomms", T_a.('Mid-Point')(in), T_a.('Has TOF')(in), count};
                    count = count + 1;
                end
            end
            
            T(count:end,:) = [];
            
            % -- Table Formatting --
            T = [T(:,2), T(:,1), T(:,3:end)];                % Sord by Date
            T = sortrows(T);
            
            T.('#') = (1:count-1)';                         % Add column to indicate mission number
            T = [T(:,end), T(:,2), T(:,1), T(:,3:end-1)];   % Return table to origonal varialbe order but with missiton number first
            
            obj.manifest = T;
            
            if dispIO, obj.Disp_Manifest; end
 
        end
        
        
        % Calculate the true vehicle speed from GPS and time stamp data
        function obj = Model_VehicleSpeed(obj)
            
            % Assume Mission Has been selected
            speed     = obj.data(1).attitude(:,5);
            timestamp = obj.data(1).timeStamp;
            lat_lon   = obj.path.GT.lat_lon;
            
            timestep = seconds(timestamp(2:end) - timestamp(1: end-1));
            timestep = [0;timestep];
            
            % filter out times where the vehicle is not moving
            out = speed == 0;
            timestep(out)  = [];
            lat_lon(out,:) = [];
            
            if isfield(obj.vehicleModel, 'speedFit')
                
                obj.filteredData.speed    = speed;
                obj.filteredData.speed(:) = obj.vehicleModel.speedFit;
                obj.data.timeStep         = timestep;
                
            else
                
                
                % Covert locations into utms for representation in meters
                [xx, yy, ~] = deg2utm(lat_lon(:,1), lat_lon(:,2));
                
                % Calculate the distance traveled between each time step
                dists = vecnorm( [xx(2:end), yy(2:end)] - [xx(1:end-1), yy(1:end-1)] , 2, 2);
                
                % Calculate the actual speed of the vehicle
                actualSpeed = dists./timestep(2:end);
                
                % Get the avarage speed and standard deviation
                mu  = mean(actualSpeed);
                sig = std(actualSpeed);
                
                speed(speed > 0) = mu;
                
                obj.vehicleModel.speed       = makedist('Normal', 'mu', 0, 'sigma', sig);
                obj.vehicleModel.units.speed = "m/s";
                obj.filteredData.speed       = speed;
                obj.data.timeStep            = timestep;
                
            end
        end
        
        
        % Calculate the altimiter Noise
        function obj = Model_Altimiter(obj, plotIO)
            
%             mapAlt = obj.orical.gtAlt;
            
            if isfield(obj.orical, 'drAlt')
                mapAlt = obj.orical.drAlt;
            else
                mapAlt = [16, 0];
            end
            
            alt = abs( obj.data.bathymetry(:,3));
            
            lat_lon = obj.data.bathymetry(:,1:2);
            
            [xx,yy] = deg2utm(lat_lon(:,1), lat_lon(:,2));
            
            dx = xx - obj.path.GT.utm(:,1);
            dy = yy - obj.path.GT.utm(:,2);
            
            obj.filteredData.alt_offset = [dx,dy];
            
%             alt = obj.data(1).attitude(:,4);
            alt_raw = alt;

            
%             if size(alt) ~= size(mapAlt) && numel(mapAlt) ~= 2
%                 warning('Altimiter profile and bathymetric profile are diferent sizes')
%             end
            
            
            out = alt > max(mapAlt(:)) | ...
                alt < 0 | ... min(mapAlt(:)) | ...
                alt == 0;
            
            sig_ = nanstd(alt(~out));
            alt(out) = 0;
            
%             alt_raw(out) = NaN;
            
            sig = sig_;
            x   = alt(1);
            for ii = 1:numel(alt)
                
                if alt(ii) == 0
                    sig = sig * 1.1;
                    continue
                end
                
                z = alt(ii);
                
                A = 1;                                  % Process Model
                B = 1;                                  % Control Model
                C = 1;                                  % Observation Matrix
                
                Q = 1;                                  % Uncertainty in the process model
                R = sig;                                % Measurment Standard Deviation
                
                u = 0;                                  % Control input
                
                
                [x, sig] = Kalman_Filter(x, sig, u, z, A, B, C, Q, R);        % Filter the measurment
                
                
                alt(ii) = x;
                
            end
            
            if isfield(obj.vehicleModel, 'altFit')                          % if the vehicle has been modeled than use that modle to correct the altimiter
                
%                 alt(alt > 0) = alt(alt > 0) + obj.vehicleModel.altFit;
                
            else
                
                obj.vehicleModel.altimeter = makedist('Normal','mu', 0,'sigma',sig_ * 3);
                obj.vehicleModel.units.altimeter = "m";
                
                
            end
            
            obj.filteredData.altitude = alt;
            obj.filteredData.altSig   = sig;
            
            
            if nargin == 2 && plotIO
                
                alt(alt == 0) = NaN;
                
                figure('name', sprintf('%s Altimiter Profiles', obj.name))
                p1 = plot(alt_raw);
                hold on
                p2 = plot(alt); % obj.orical.gtAlt);
                hold off
                
                xlabel('Time Step')
                ylabel('Altitude [m]')
                title( sprintf('%s Altimiter Profile', obj.name) )
                legend([p1,p2], {'Raw Altimeter Data', 'Filtered Altimiter Data'}) % GT Bathy Profile'} )  %, 'DR Bathy Profile'
            end
            
        end
        
        
        % Create model of the compass
        function obj = Model_Compass(obj, modelType, plotIO)
            
            compass = obj.data(1).attitude(:,3);                            % Comapss heading
            
            if isfield(obj.vehicleModel, 'compassFit')
                
                f = obj.vehicleModel.compassFit;
                obj.filteredData.heading = feval(f,compass);
                
            else
                
                gps = obj.path.GT.utm;                                         % Get Vehicle Location
%                 gps = zeros(size(gps0));
%                 
%                 % Smooth out the GPS data to determin the true vehicle heading
%                 n = size(gps0,1);
%                 sw = SlidingWindow(20, n);
%                 
%                 for ii = 1:n
%                     
%                     data_ = gps0(sw.D1(ii),:);
%                     gps1(ii,:) = mean(data_);
%                     
%                 end
                
                % Determin the true vehicl;e heading
                dxdy = zeros(size(gps));
                dxdy(1:end-1,:) = gps(2:end, :) - gps(1:end-1,:);
                dxdy(end,:) = dxdy(end-1,:);
                
                oo = dxdy(:,1) == 0 & dxdy(:,2) == 0;
                
                angle = NaN(size(compass));
                angle(~oo) = atan2(dxdy(~oo,2), dxdy(~oo,1)) * 180/pi;          % Deturmine true heading
                
                neg = angle < 0;                                                % Keep angles between 0 and 360
                angle(neg) = angle(neg) + 360;
                
                angle = CompassAngle(angle, 'deg', 'deg');                      % Conver to compass angle
                
                
                % Get Compass Error and compensate for rolle over
                err = angle - compass;                                          % Compass Error
                
                for ii = find( abs(err) > 180)'
                    
                    if angle(ii) > 180
                        angle(ii) = angle(ii) - 360;
                    else
                        angle(ii) = angle(ii) + 360;
                        
                    end
                end
                
                err = angle - compass;                                          % Compass Error
                
                angle(isnan(angle)) = compass(isnan(angle));
                
                %             mu = nanmean(err);
                sig = nanstd(err);
                
                obj.vehicleModel.compass = makedist('Normal','mu', 0,'sigma',sig);
                obj.vehicleModel.units.compass = 'deg';
                
                c = compass(~isnan(angle));
                a = angle(~isnan(angle));
                
                f = fit( c, a, modelType);
                
                y = feval(f,compass);
                
                obj.filteredData.heading = y;
                
                
                if nargin == 3 && plotIO
                    figure('name', sprintf("Compass Calibration for %s",obj.name), 'numbertitle', 'off');
                    p1 = plot(angle, 'g');
                    hold on
                    p2 = plot(compass, 'b');
                    p3 = plot(err, 'r');
                    p4 = plot(y, 'm');
                    hold off
                    
                    xlabel('Time Steps', 'FontSize', 14)
                    ylabel('Compass Heading [deg]', 'FontSize', 14)
                    title(sprintf("Compass Calibration for %s",obj.name), 'FontSize', 14)
                    legend([p1,p2,p3,p4], {'gps','compass','errer','fit'}, 'FontSize', 12)
                    drawnow
                end
                
            end
        end
        
        
        
        % Select data from an indivual mission
        function obj = Select_Mission(obj, idx)
            
            % Convert time string into datetime
%             if nargin == 3, missionTime = datetime(timeStr, 'Format','dd-MMM-yyyy HH:mm:ss'); end
            
            % Reset data objects
            obj.data   = obj.rawData;
            obj.acomms = obj.rawAcomms;
            obj.course = obj.rawCourse;
            
            % Get Rid of old models
            obj.filteredData = [];
%             obj.vehicleModel = [];
            
            mission    = obj.manifest.Mission(idx);
            start_time = obj.manifest.('Start Time')(idx);
            end_time   = obj.manifest.('End Time')(idx);
            
            % Iterate through the data objects to find the desired mission
            for type = {'rawData', 'rawAcomms', 'rawCourse'; ...
                        'data',    'acomms',    'course'}
                
                for bb = size(obj.(type{1}),1): -1:1
                    
                    % --- See if the data set contains the desired mission
                    if strcmp(type{1}, 'rawCourse')
                        in = strcmp( obj.(type{1})(bb).mission, mission);
                    else
                        times = obj.(type{1})(bb).gpsDate;
                        in = ismember( obj.(type{1})(bb).mission, mission) ... 
                             & start_time <= times & times <= end_time;
                    end
                    
                    % --- Mission or time is not present --> remove entry from data sets
                    if ~any(in)
                        obj.(type{2})(bb) = [];
                        continue
                    end
                    
                    
                    % --- Mission (and start time) are present, --> only retain the data lines corresponding to the mission
                    for mems = fields(obj.(type{2})(bb))'
                        
                        if strcmp( mems, 'header') || strcmp( mems, 'log') || strcmp( mems, 'timeErr')  % Skip the headers and log entrances
                            continue
                        elseif strcmp( mems, 'mission') || strcmp( mems, 'waypoints')                   % Copy over the entire mission name and waypoints
                            obj.(type{2})(bb).(mems{1}) = obj.(type{1})(bb).(mems{1});
                        else                                                                            % Choose individual data that are associated witht he mission
                            data_ = obj.(type{1})(bb).(mems{1});
                            obj.(type{2})(bb).(mems{1}) = data_(in,:);
                            
                        end
                        
                    end
               
                end
            end
            
            
            obj.filteredData = [];                                          % Get Rid of filtered Data
        
        end
        
        
        function [pf, path, name_, r, speed, speedNoise, compassNoise] = VehicleInfo(obj, name_)
            
            if nargin == 1
                name_ = obj.name;
            end
            
            pf = obj.tbn;
            r  = 3;
            
            speed = obj.filteredData.speed(1);
            
            speedNoise   = obj.vehicleModel.speed;
            compassNoise = obj.vehicleModel.compass;
            
            path = fliplr(obj.course.waypoints);
            
            trueStart = obj.path.GT.lat_lon(1,:);
            
            dxdy = path(1,:) - trueStart;
            
            path = path - dxdy;
            
            try
            [xx, yy] = deg2utm( path(:,1), path(:,2) );
            catch E
                disp(E)
            end
            
            path = [xx,yy];
            
            
            % Pack the information into a cell array so that we can nest this function call inside another function call
            if nargout == 1
            
                pf = {pf, path, name_, r, speed, speedNoise, compassNoise};
                
            end
            
            
        end
        
        
        function obj = VehicleModeling(obj, modelType, bathy, skip, dispIO)
            
            a_     = zeros(10000,1);
            c_     = zeros(10000,1);
            err_   = zeros(10000,1);
            speed_ = zeros(10000,1);
            altEr_ = zeros(10000,1);
            
            a = 1;
            for ii = 1:size(obj.manifest,1)
   
                if nargin >= 4 && any(ii == skip), continue, end
                
                
                obj = obj.Select_Mission(ii);                               % Get Mission
                obj = obj.Get_VehiclePaths;                                 % Get smoothed GPS path
                obj = obj.Model_Altimiter;                                  % Model altimiter
                
                % --- Model Vehicle Speed ---------------------------------
                speed     = obj.data(1).attitude(:,5);
                compass   = obj.data(1).attitude(:,3);                        % Comapss heading
                timestamp = obj.data(1).timeStamp;
                alt       = obj.filteredData.altitude;
                utm       = obj.path.GT.utm;
                
                
                % filter out times where the vehicle is not moving
                out = speed == 0;
                timestamp(out) = [];
                compass(out)   = [];
                utm(out,:)     = [];
                alt(out,:)     = [];
                
                % Calculate the actual speed of the vehicle
                timestep = seconds( timestamp(2:end) - timestamp(1: end-1) );
                dists    = vecnorm( utm(2:end,:) - utm(1:end-1,:) , 2, 2 );
                
                b = a+ numel(dists) -1;
                
                speed_(a:b) = dists./timestep;
                
                
                % --- Model Vehicle Compass -------------------------------
                % Determin the true vehicle heading
                dxdy = utm(2:end,:) - utm(1:end-1,:);
                compass(end)   = [];
                
%                 if any(dxdy == 0), warning("0s in dx -dy"), end

                
                
                angle = atan2(dxdy(:,2), dxdy(:,1)) * 180/pi;               % Deturmine true heading
                
                neg = angle < 0;                                            % Keep angles between 0 and 360
                angle(neg) = angle(neg) + 360;
                
                angle = CompassAngle(angle, 'deg', 'deg');                  % Conver to compass angle
                
                err = angle - compass;                                      % Compass Error
                
                for jj = find( abs(err) > 180)'                             % Correct for rolle over
                    
                    if angle(jj) > 180
                        angle(jj) = angle(jj) - 360;
                    else
                        angle(jj) = angle(jj) + 360;
                    end
                end
                
                a_(a:b)   = angle(~isnan(angle));
                c_(a:b)   = compass(~isnan(angle));
                err_(a:b) = angle - compass;                                % Compass Error
                
                
                % --- Modle Altimiter -------------------------------------
                gtAlt = bathy.Depth(utm);
                
                gtAlt(end) = [];
                alt(end)   = [];
                
                alt(alt == 0) = NaN;
                
                altEr_(a:b) = gtAlt - alt;
                
                
                
                a = b+1;
                
            end
            
            
            % Get Rid of extra entries
            a_(a:end)     = [];
            c_(a:end)     = [];
            err_(a:end)   = [];
            speed_(a:end) = [];
            
            % Get rid of zerso
            a_(a_ == 0) = [];
            c_(c_ == 0) = [];
            err_(err_ == 0)   = [];
            speed_(speed_ == 0) = [];
            
            
            % Speed Model
            mu  = mean(speed_);
            sig = std(speed_);
           
            stat.speed.mu = mu;
            stat.speed.sig = sig;
            
            obj.vehicleModel.speed       = makedist('Normal', 'mu', 0, 'sigma', sig);
            obj.vehicleModel.units.speed = "m/s";
            
            obj.vehicleModel.speedFit = mu;
            
            
            
            % Compass model
            mu = nanmean(err);
            sig = nanstd(err_);
            
            stat.compass.mu = mu;
            stat.compass.sig = sig;
            
            obj.vehicleModel.compass = makedist('Normal','mu', 0,'sigma',sig);
            obj.vehicleModel.units.compass = 'deg';
            
            obj.vehicleModel.compassFit     = fit( c_, a_, modelType);
            obj.vehicleModel.compassFit_inv = fit( a_, c_, modelType);
            
            
            % Altimiter Model
            mu  = nanmean( altEr_ );
            sig = nanstd( altEr_ );
            
            stat.altimeter.mu = mu;
            stat.altimeter.sig = sig;
            
            obj.vehicleModel.altimeter = makedist('Normal','mu', 0,'sigma',sig);
            obj.vehicleModel.units.altimeter = "m";
            
            obj.vehicleModel.altFit = mu;
            
            
            if nargin == 5 && dispIO
                
                fprintf('\n%s Vehicle Stats:\n', obj.name)
                
                for ii = {'Altimeter', 'Compass', 'Speed';
                          'altimeter', 'compass', 'speed';
                          '[m]',       '[deg]',   '[m/s]'}
                    
                    mu  = stat.(ii{2}).mu;
                    sig = stat.(ii{2}).sig;
                    
                    fprintf("%s: %.3f %s %.3f  %s\n",ii{1}, mu, char(177), sig, ii{3})
                    
                end
            end
            
        end
         
        
        % _____ Plotting __________________________________________________________________________
        function figOut = Plot_Paths(obj, fig_, paths, c)
            
            
            if nargin > 1 && isa(fig_, 'Geotiff')
                fig = figure('name', sprintf('Riptide %s Paths', obj.name), 'numbertitle', 'off');
                fig_.Show
                
            elseif nargin > 1 && isgraphics(fig_, 'figure')
                fig = figure(fig_);
            else
                fig = figure('name', sprintf('Riptide %s Paths', obj.name), 'numbertitle', 'off');
            end
            
            
            if nargin < 3, mem = fields(obj.path(1));
            else,          mem = paths';
            end
            
            if nargin < 4, c = 'grbcm'; end
            
            hold on
            
            
            if ismember('EKF', mem)
            
                xx = obj.path.EKF.utm(:,1);
                yy = obj.path.EKF.utm(:,2);
                zone = obj.path.GT.utmzone(1,:);
                
                
                [lat, lon] = utm2deg(xx,yy, repmat(zone, size(xx,1),1) );
                
                obj.path.EKF.lat_lon = [lat, lon];
                
            end
            
            for ii = numel(mem): -1: 1
                
                lat_lon = obj.path(1).(mem{ii}).lat_lon;
                
                p(ii) = plot(lat_lon(:,2),   lat_lon(:,1),   c(ii) );       % Plot Path
                p1    = plot(lat_lon(1,2),   lat_lon(1,1),   ['d',c(ii)] ); % Plot Path Begining
                p2    = plot(lat_lon(end,2), lat_lon(end,1), ['*',c(ii)] ); % Plot Path End
                
            end
            
            hold off
            
            ax = gca;
            ax.FontSize = 12;       % Set font size first or else you will loose the rest of the formatting
            
            title(sprintf("Zig Zag Mission 2: Riptide %s's Paths", obj.name), 'FontSize', 20)
            xlabel('Longitude', 'FontSize', 16)
            ylabel('Latitude',  'FontSize', 16)
            legend([p, p1, p2] ,[mem; {'Start'; 'End'}], 'Interpreter', 'non', 'FontSize', 14)
            
            drawnow
            
            if nargout == 1
                figOut = fig;
            end
            
        end
          
        
        function figOut = Plot_Course(obj, fig_)
            
            if nargin > 1 && isa(fig_, 'Geotiff')
                fig = figure('name', sprintf('Riptide %s Paths', obj.name), 'numbertitle', 'off');
                fig_.Show;
                
            elseif nargin > 1 && isgraphics(fig_, 'figure')
                fig = figure(fig_);
            else
                fig = figure('name', sprintf('Riptide %s Paths', obj.name), 'numbertitle', 'off');
            end
            
            waypoints = obj.course.waypoints;
           
            hold on
            
            plot(waypoints(:,1), waypoints(:,2), '*-c')
            plot(waypoints(1,1), waypoints(1,2), 'dc')
            
            
            hold off
            
            ax = gca;
            ax.FontSize = 12;       % Set font size first or else you will loose the rest of the formatting

            title(sprintf("Riptide %s's Course", obj.name), 'FontSize', 20)
            xlabel('Longitude', 'FontSize', 16)
            ylabel('Latitude',  'FontSize', 16)
            
            drawnow
            
            if nargout == 1
                figOut = fig;
            end
            
            
        end
        
        
        function figOut = Plot_Acomms(obj, fig_)
            
            if nargin > 1 && isa(fig_, 'Geotiff')
                fig = figure('name', sprintf('Riptide %s A-comms', obj.name), 'numbertitle', 'off');
                geoshow(fig_.Image, fig_.Data)
                
            elseif nargin > 1 && isgraphics(fig_, 'figure')
                fig = figure(fig_);
            else
                fig = figure('name', sprintf('Riptide %s A-comms', obj.name), 'numbertitle', 'off');
            end
            
            hold on
            plot(obj.acomms.vehicle(:,2), obj.acomms.vehicle(:,1), '*')
            hold off
            
            ax = gca;
            ax.FontSize = 12;       % Set font size first or else you will loose the rest of the formatting

            title(sprintf("Riptide %s's Acoustric Communications", obj.name), 'FontSize', 20)
            xlabel('Longitude', 'FontSize', 16)
            ylabel('Latitude',  'FontSize', 16)
 
            
            drawnow
            
            if nargout == 1
                figOut = fig;
            end
            
        end
        
        
        function figOut = Plot_AltimeterProfile(obj, fig)
            
            times = obj.data.timeStamp;
            
            if isfield(obj.filteredData, 'altitude')
                alt = obj.filteredData.altitude;
            else
                alt = obj.data(1).attitude(:,4);
            end
            
%             if nargin == 2
%                 figure(fig)
%                 hold on
%             else
%                 fig = figure('name',sprintf('%s Altimeter Measurments',obj.name), 'numbertitle', 'off');
%             end
%             
%             if numel(times) == numel(alt)
%                 plot(times, alt)
%             else
%                 warning("Plot_Altitude --> Times and Alt are different sizes")
%                 plot(alt)
%             end
%             
%             title(sprintf('%s: Atitude vs Time', obj.name), 'FontSize', 14)
%             xlabel('Date-Time', 'FontSize', 14)
%             ylabel('Atitude [m]', 'FontSize', 14)
%             legend('FontSize', 12)
%             drawnow
            

                alt(alt == 0) = NaN;
                
                fig = figure('name', sprintf('%s Altimiter Profiles', obj.name));
                p1 = plot(alt, '.');
                hold on
                p2 = plot(obj.orical.gtAlt);
                p3 = plot(obj.orical.drAlt);
                hold off

                ax = gca;
                ax.FontSize = 12;       % Set font size first or else you will loose the rest of the formatting
                
                xlabel('Time Step')
                ylabel('Altitude [m]')
                title( sprintf('%s Altimiter Profile', obj.name) )
                legend([p1,p2, p3], {'Filtered Alt', 'GT Bathy Profile', 'DR Bathy Profile'} )
                drawnow
                
                
            if nargout == 1
                figOut = fig;
            end
            
        end
        
        
        function figOut = Plot_Altimeter(obj, fig_)
            
            
            if nargin > 1 && isa(fig_, 'Geotiff')
                fig = figure('name', sprintf("Riptide %s's Altimiter Data", obj.name), 'numbertitle', 'off');
                fig_.Show
                
            elseif nargin > 1 && isgraphics(fig_, 'figure')
                fig = figure(fig_);
            else
                fig = figure('name', sprintf("Riptide %s's Altimiter Data", obj.name), 'numbertitle', 'off');
            end
            
            
            if isfield(obj.filteredData, 'altitude')
                alt = obj.filteredData.altitude;
            else
                alt = obj.data(1).attitude(:,4);
            end
            
            oo = alt == 0;
            
            lat_lon = obj.path.GT.lat_lon;
            
            hold on
            
            p1 = plot(lat_lon(oo,2), lat_lon(oo,1), 'r.');        % Plot 0s
            p2 = plot(lat_lon(~oo,2), lat_lon(~oo,1), 'g.');      % Plot ~0s
            
            hold off

            ax = gca;
            ax.FontSize = 12;       % Set font size first or else you will loose the rest of the formatting
            
            title(sprintf("Riptide %s's Altimiter Data", obj.name), 'FontSize', 14)
            xlabel('Easitng', 'FontSize', 14)
            ylabel('Northing', 'FontSize', 14)
            legend([p1, p2] ,{'Zeros'; 'Valid Data'}, 'Interpreter', 'non', 'FontSize', 12)
            
            drawnow
            
            if nargout == 1
                figOut = fig;
            end
            
        end
        
        
        function figOut = Show_Manifest(obj)
            
            obj.Disp_Manifest;
            
            fig = figure('name','Riptide Manifest','numbertitle','off');
            
            % Get the table in string form.
            TString = evalc('disp(obj.manifest)');
            % Use TeX Markup for bold formatting and underscores.
            TString = strrep(TString,'<strong>','\bf');
            TString = strrep(TString,'</strong>','\rm');
            TString = strrep(TString,'_','\_');
            % Get a fixed-width font.
            FixedWidth = get(0,'FixedWidthFontName');
            % Output the table using the annotation command.
            annotation(gcf,'Textbox','String',TString,'Interpreter','Tex',...
                'FontName',FixedWidth,'Units','Normalized','Position',[0 0 1 1]);
            
            
            %             uit = uitable(fig, 'Data',T);
            %             uit.Units =  'Normalized';
            %             uit.Position = [20, 20, 100, 100];
            
            if nargout == 1
                figOut = fig;
            end
            
        end
        
        
        function figOut = Plot_TimeStamps(obj)
            
            times = obj.data.timeStamp;
            
            fig = figure('name', sprintf('%s Mission Timestamps', obj.name), 'numbertitle', 'off');
            plot(times)
            
            title(sprintf('%s Mission Timestamps', obj.name), 'FontSize', 14)
            xlabel('Time Step', 'FontSize', 14)
            ylabel('Date-Time', 'FontSize', 14)
            drawnow
            
            if nargout == 1
                figOut = fig;
            end
            
        end
        
        
    end
    
    
    methods (Access = private)
        
        % Calculate Speed and Heading Error
        function[errTheta, errSpeed] = GetVelocityErrors(gpsX, gpsY, speed, heading, timestamp)
            
            % Get Path from Speed and heading
            timeStep = seconds(timestamp(2:end) - timestamp(1: end-1));
            
            deltaX_dr = [0; speed(1:end-1) .* timeStep .* sin( heading(1:end-1) )];
            deltaY_dr = [0; speed(1:end-1) .* timeStep .* cos( heading(1:end-1) )];
            
            deltaX_gps = [0; gpsX(2:end) - gpsX(1:end-1)];
            deltaY_gps = [0; gpsY(2:end) - gpsY(1:end-1)];
            
            % Error between deadreckoning and GPS
            errX = deltaX_dr - deltaX_gps;
            errY = deltaY_dr - deltaY_gps;
            
            
            % Deturmin uncertainties ---------------------------------------------------------------------
            sigtrig = sin(heading).^2 - cos(heading).^2;
            
            sigX_dX = errX ./ deltaX_dr;
            sigY_dY = errY ./ deltaY_dr;
            
            errTheta = abs(heading) .* sqrt( abs((sigY_dY .* sigY_dY  - sigX_dX .* sigX_dX) ./ sigtrig ));
            errSpeed = abs(speed)   .* sqrt( abs((tan(heading ).^2 .*  sigX_dX .^2 - sigY_dY.^2) ./ ( tan(heading ).^2 - 1)));
            
            errTheta(isinf(errTheta)) = [];
            errSpeed(isinf(errSpeed)) = [];
            
            errTheta = nanmean(errTheta);
            errSpeed = nanmean(errSpeed);
            
            fprintf("\n\nTheta Err: %f\tSpeed Err: %f\n\n",errTheta, errSpeed)
            
        end
        
        
        % Plotting Template
        function figOut = Plot_Template(obj, fig_)
            
            
            if nargin > 1 && isa(fig_, 'Geotiff')
                fig = figure('name', sprintf("Riptide %s's Altimiter Data", obj.name), 'numbertitle', 'off');
                fig_.Show
                
            elseif nargin > 1 && isgraphics(fig_, 'figure')
                fig = figure(fig_);
            else
                fig = figure('name', sprintf("Riptide %s's Altimiter Data", obj.name), 'numbertitle', 'off');
            end
            
            
            if isfield(obj.filteredData, 'altitude')
                alt = obj.filteredData.altitude;
            else
                alt = obj.data(1).attitude(:,4);
            end
            
            oo = alt == 0;
            
            lat_lon = obj.path.GT.lat_lon;
            
            
            p1 = plot(lat_lon(oo,2), lat_lon(oo,1), 'r.');        % Plot 0s
            p2 = plot(lat_lon(~oo,2), lat_lon(~oo,1), 'g.');      % Plot ~0s
            
            hold off
            
            title(sprintf("Riptide %s's Altimiter Data", obj.name), 'FontSize', 14)
            xlabel('Easitng', 'FontSize', 14)
            ylabel('Northing', 'FontSize', 14)
            legend([p1, p2] ,{'Zeros'; 'Valid Data'}, 'Interpreter', 'non', 'FontSize', 12)
            
            drawnow
            
            if nargout == 1
                figOut = fig;
            end
            
        end
        
        
    end
end




function Data = GetData(file, fun, fileType)

if isfile(file), Data = fun(file);
    
elseif isfolder(file), Data = BatchDir(file,fileType,fun);
    
else, disp("RT_Data --> Invalid File Path")
    
end

end



function course = GetWaypoints(file)

[~, fileName, ext] = fileparts(file);

fprintf(' * Importing Riptide Mission: %s\n', [fileName,ext]);

try
    fileID = fopen(file,'r');
    
    line = "";
    
    while ~strncmp(line, 'Mission:', 8), line = fgetl(fileID); end
    
    mission = line;
    
    while ~strcmp(line, 'Waypoints:'), line = fgetl(fileID); end
    
    dataArray = textscan(fileID, '%f%f', 'Delimiter', ',', 'ReturnOnError', false);
    
    fclose(fileID);
    
    
    course.waypoints = [dataArray{1}, dataArray{2}];
    
    mission = strsplit(mission,':');
    course.mission = strtrim(mission{2});
    
    
catch exception                 % File Cannot be read
    warning off backtrace
    warning('Error Reading waypoint file \n\t\t%s \n\t\t%s \n\t\tError in: %s, Line: %d \n\n\t\tSuspect bad log file, skipping it \n\n',exception.identifier,exception.message,exception.stack(1).name,exception.stack(1).line)
    warning on backtrace
    
    
end

end


function T = missionTabe(data, skipIdel)

is_Acomms = isfield(data, 'tof');

if is_Acomms
    sz = [100,5];
    varTypes = {'string', 'datetime',  'datetime', 'datetime',  'logical'};
    varNames = {'Mission','Start Time','End Time', 'Mid-Point', 'Has TOF'};
    
    tof = cat(1, data.tof);
    
else
    sz = [100,4];
    varTypes = {'string', 'datetime',  'datetime', 'datetime'};
    varNames = {'Mission','Start Time','End Time', 'Mid-Point'};
    
end

T = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);

time_data = cat(1,data.gpsDate);
mis_all   = cat(1,data.mission);
mis_unq = unique(mis_all);

% sourt the data chronologicaly to avoid funny-ness from overlappint mission types
[time_data, I] = sort(time_data);

mis_all = mis_all(I);

t_idx = 1;
% Find the number of times that a mission was run and the start and end times
for mis = mis_unq'
    
    if skipIdel && strcmp(mis, 'idel'),    continue, end
    
    is_in = false;
    count = 0;
    time_start = NaT(100,1);
    time_end   = NaT(100,1);
    
    indx = NaN(100,2);
    
    iter  = 1;
    
    for in = strcmp(mis_all, mis)'
        
        if in && ~is_in
            count = count +1;
            time_start(count) = time_data(iter);
            indx(count,1) = iter;
            
        elseif ~in && is_in
            time_end(count) = time_data(iter-1);
            indx(count,2) = iter -1;
            
        end
        is_in = in;
        iter = iter+1;
        
    end
    
    if isnat(time_end(count))
        time_end(count) = time_data(iter-1);
        indx(count,2) = iter-1;
    end             % End time of the last mission
    
    
    if is_Acomms
        for ii = 1:count
            
            tof_ = isnan( tof(indx(ii,1) : indx(ii,2), :) );
            tof_ = ~all(tof_(:));
            
            T(t_idx,:) = {mis{1} ,time_start(ii), time_end(ii), time_start(ii) + (time_end(ii) - time_start(ii))/2, tof_};
            t_idx = t_idx+1;
        end
        
    else
        
        
        for ii = 1:count
            T(t_idx,:) = {mis{1} ,time_start(ii), time_end(ii), time_start(ii) + (time_end(ii) - time_start(ii))/2};
            t_idx = t_idx+1;
        end
        
    end
end

T(t_idx:end,:) = [];

end



