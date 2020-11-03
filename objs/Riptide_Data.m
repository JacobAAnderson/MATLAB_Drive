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
        path
        tbn
        rawAcomms
        rawCourse
        rawData
        vehicleModel
        pAcommsHandler
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
        
        
        % Read in mission files to get the waypoints
        function obj = Add_Waypoints(obj, file)
            
            fun = @(file_) GetWaypoints(file_);
            obj.rawCourse = GetData(file, fun, 'txt');
            obj.course = obj.rawCourse;
            
        end
        
        
        % Provide the vehicle with a particle filter for Terrain Based Navigation annalisys
        function obj = Add_TBN(obj, map)
            
            obj.tbn = DEC_TBN(obj.name, 1000 );                              % Instanciate Paricle Filter
            
            obj.tbn = obj.tbn.Add_Map(map);                                   % Provide Bathymetry map
            obj.tbn = obj.tbn.SetLocation(obj.path.GT.utm(1,:), [3,0;0,3]);   % Tell it where the vehicle is starting from
            
 
            % --- Set vehicle transiton model as a probability distribuition ( "type of distribuition' , mu, sigma ) ---
            for prop = {'acoms',  'altimeter', 'compass' 'speed';           % --> Type of noise
                        'Normal', 'Normal',    'Normal', 'Normal';          % --> Default distribion type
                         0,        0,           0,        0;                % --> Default mu
                         0.5,      0.5,         10,       1.0}              % --> Default sigma
                     
                     
                     if  isfield( obj.vehicleModel, prop{1} )
                         mu     = obj.vehicleModel.(prop{1}).mu;
                         sig    = obj.vehicleModel.(prop{1}).sigma;
                         name_  = obj.vehicleModel.(prop{1}).DistributionName;
                         obj.tbn = obj.tbn.SetNoise(prop{1}, name_, mu, sig);
                         
                     else
                         obj.tbn = obj.tbn.SetNoise(prop{1}, prop{2}, prop{3}, prop{4});
                     end
                     
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
      
        
        % Get Vehicle Path in UTM an degrees lat-lon
        function obj = Get_VehiclePaths(obj)
            
            obj.path = [];
            
            speed     = obj.data(1).attitude(:,5);
            heading   = obj.data(1).attitude(:,3);
            timestamp = obj.data(1).timeStamp;
            
            timeStep_ = [0; seconds(timestamp(2:end) - timestamp(1: end-1))];
            
            deltaX_dr = speed .* timeStep_ .* sind( heading );
            deltaY_dr = speed .* timeStep_ .* cosd( heading );
            
            x = cumsum(deltaX_dr);
            y = cumsum(deltaY_dr);
            
            [xx, yy, utmzone] = deg2utm( obj.data(1).vehicle(:,1), obj.data(1).vehicle(:,2) );
            
            [dr_lat, dr_lon] = utm2deg(x + xx(1), y+yy(1), utmzone);
            
            obj.path.GT.lat_lon = obj.data(1).vehicle(:,1:2);
            obj.path.GT.utm     = [xx,yy];
            obj.path.GT.utmzone = utmzone;
            
            obj.path.DR.utm(:,1) = x + xx(1);
            obj.path.DR.utm(:,2) = y + yy(1);
            obj.path.DR.lat_lon  = [dr_lat, dr_lon];
            
            obj.data.timeStep = timeStep_;
            
            if ~ ( isfield( obj.vehicleModel, 'speed') || isfield( obj.filteredData, 'heading') ), return, end
            
            if isfield( obj.vehicleModel, 'speed'),   speed   = obj.filteredData.speed;   end
            if isfield( obj.filteredData, 'heading'), heading = obj.filteredData.heading; end
            
            deltaX_dr = speed .* timeStep_ .* sind( heading );
            deltaY_dr = speed .* timeStep_ .* cosd( heading );
            
            x = cumsum(deltaX_dr);
            y = cumsum(deltaY_dr);
            
            [xx, yy, utmzone] = deg2utm( obj.data(1).vehicle(:,1), obj.data(1).vehicle(:,2) );
            
            [dr_lat, dr_lon] = utm2deg(x + xx(1), y+yy(1), utmzone);
            
            obj.path.DR_corrected.utm(:,1) = x + xx(1);
            obj.path.DR_corrected.utm(:,2) = y + yy(1);
            obj.path.DR_corrected.lat_lon  = [dr_lat, dr_lon];
            
            
            
        end
        
        
        % Retrive the path recorded by the particle filter
        function obj = Get_TBNpath(obj)
            
            utm = obj.tbn.path(1:obj.tbn.path_idx-1,1:2);
            
            zone = obj.tbn.Map.UTM_Zone;
            
            [lat, lon]  = utm2deg( utm(:,1), utm(:,2), repmat(zone, size(utm,1),1 )) ;
            
            obj.path.tbn.utm = utm;
            obj.path.tbn.lat_lon = [lat, lon];
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
            lat_lon   = obj.data(1).vehicle(:,1:2);
            timestamp = obj.data(1).timeStamp;
                
            timestep = seconds(timestamp(2:end) - timestamp(1: end-1));
            timestep = [0;timestep];
            
            % filter out times where the vehicle is not moving
            out = speed == 0;
            timestep(out)  = [];
            lat_lon(out,:) = [];
           
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
        
        
        % Calculate the altimiter Noise
        function obj = Model_Altimiter(obj, threshold)
            
            alt = obj.data(1).attitude(:,4);
            
%             mu  = mean(alt);
%             sig = std(alt);
            
            % --- Filter Altimiter Data ----
            if nargin == 2
                alt(alt > threshold) = 0;
            end
            
            alt(alt < 0) = 0;
            alt1 = alt;
             
            A = 1;              % Process Model
            B = 1;              % Control Model
            C = 1;              % Observation Matrix
            
            Q = 1;            % Uncertainty in the process model
            R =std(alt);        % Measurment Standard Deviation
            
            u = 0;              % Control input
            
            sw = SlidingWindow(20, numel(alt));                             % Slididng window ojbect
            
            for c = 1:2
                
                x = alt(1);         % State
                sig = 0.75;         % Error in the state
                
                for ii = 2:numel(alt)
                    
                    rows = sw.D1(ii);
                    window = alt(rows);
                    oo = window == 0;
                    
                    if sum(oo) > 12, alt1(ii) = 0; continue %, end
                        
                    elseif alt(ii) == 0, z = mean(window(window ~= 0));
                    else,            z = alt(ii);
                    end
                    
                    %                 z = alt(ii);
                    
                    [x, sig] = Kalman_Filter(x, sig, u, z, A, B, C, Q, R);
                    
                    alt1(ii) = x;
                end
                
                alt = alt1;
            end
            
            obj.vehicleModel.altimeter = makedist('Normal','mu', 0,'sigma',0.75);
            obj.vehicleModel.units.altimeter = "m";
            obj.filteredData.altitude = alt1;
            
        end
        
        
        % Create model of the compass
        function obj = Model_Compass(obj, modelType, plotIO)
            
            compass = obj.data(1).attitude(:,3);                            % Comapss heading
            
            gps0 = obj.path.GT.utm;                                         % Get Vehicle Location
            gps1 = zeros(size(gps0));
            
            % Smooth out the GPS data to determin the true vehicle heading
            n = size(gps0,1);
            sw = SlidingWindow(20, n);
            
            for ii = 1:n
                
                data_ = gps0(sw.D1(ii),:);
                gps1(ii,:) = mean(data_);
                
            end
            
            % Determin the true vehicl;e heading
            dxdy = zeros(size(gps1));
            dxdy(1:end-1,:) = gps1(2:end, :) - gps1(1:end-1,:);
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
            
            legend([p1,p2,p3,p4], {'gps','compass','errer','fit'})
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
            obj.vehicleModel = [];
            
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
        
        
        
        % _____ Plotting __________________________________________________________________________
        function fig = Plot_Paths(obj, fig_)
            
            
            if nargin > 1 && isa(fig_, 'Geotiff')
                fig = figure('name', sprintf('Riptide %s Pathes', obj.name), 'numbertitle', 'off');
                geoshow(fig_.Image, fig_.Data)
                
            elseif nargin > 1 && isgraphics(fig_, 'figure')
                fig = figure(fig_);
            else
                fig = figure('name', sprintf('Riptide %s Pathes', obj.name), 'numbertitle', 'off');
            end
            
            mem = fields(obj.path(1));
            
            hold on
            c = 'grbcm';
            
            for ii = numel(mem): -1: 1
                
                lat_lon = obj.path(1).(mem{ii}).lat_lon;
                
                p(ii) = plot(lat_lon(:,2),   lat_lon(:,1),   c(ii) );       % Plot Path
                p1    = plot(lat_lon(1,2),   lat_lon(1,1),   ['d',c(ii)] ); % Plot Path Begining
                p2    = plot(lat_lon(end,2), lat_lon(end,1), ['*',c(ii)] ); % Plot Path End
                
            end
            
            hold off
            
            title(sprintf("Riptide %s's Pathes", obj.name))
            xlabel('Easitng')
            ylabel('Northing')
            legend([p, p1, p2] ,[mem; {'Start'; 'End'}], 'Interpreter', 'non')
            
            drawnow
        end
        
        
        function fig = Plot_Path(obj, fig_, path, c)
            
            if nargin > 1 && isa(fig_, 'Geotiff')
                fig = figure('name', sprintf('Riptide %s Pathes', obj.name), 'numbertitle', 'off');
                geoshow(fig_.Image, fig_.Data)
                
            elseif nargin > 1 && isgraphics(fig_, 'figure')
                fig = figure(fig_);
            else
                fig = figure('name', sprintf('Riptide %s Pathes', obj.name), 'numbertitle', 'off');
            end
            
            if nargin < 4, c = 'g'; end
            
            ii = 1;
            
            figure(fig);
            hold on
            
            plot(obj.path(ii).(path).lat_lon(:,2),   obj.path(ii).(path).lat_lon(:,1),   c);     % Plot Pathe
            plot(obj.path(ii).(path).lat_lon(end,2), obj.path(ii).(path).lat_lon(end,1), ['*',c]);
            plot(obj.path(ii).(path).lat_lon(end,2), obj.path(ii).(path).lat_lon(end,1), ['o',c]);
            
            hold off
            drawnow
        end
        
        
        function fig = Plot_Course(obj, fig_)
            
            if nargin > 1 && isa(fig_, 'Geotiff')
                fig = figure('name', sprintf('Riptide %s Pathes', obj.name), 'numbertitle', 'off');
                geoshow(fig_.Image, fig_.Data)
                
            elseif nargin > 1 && isgraphics(fig_, 'figure')
                fig = figure(fig_);
            else
                fig = figure('name', sprintf('Riptide %s Pathes', obj.name), 'numbertitle', 'off');
            end
            
            waypoints = obj.course.waypoints;
            
            str = fig.CurrentAxes.Legend.String;
            hold on
            
            plot(waypoints(:,1), waypoints(:,2), '*-c')
            plot(waypoints(1,1), waypoints(1,2), 'dc')
            
            fig.CurrentAxes.Legend.String = str;
            
            hold off
            drawnow
        end
        
        
        function fig = Plot_Acomms(obj, fig_)
            
            if nargin > 1 && isa(fig_, 'Geotiff')
                fig = figure('name', sprintf('Riptide %s Pathes', obj.name), 'numbertitle', 'off');
                geoshow(fig_.Image, fig_.Data)
                
            elseif nargin > 1 && isgraphics(fig_, 'figure')
                fig = figure(fig_);
            else
                fig = figure('name', sprintf('Riptide %s Pathes', obj.name), 'numbertitle', 'off');
            end
            
            hold on
            plot(obj.acomms.vehicle(:,2), obj.acomms.vehicle(:,1), '*')
            hold off
            drawnow
        end
        
        
        function fig = Plot_Altitude(obj, fig)
            
            times = obj.data.timeStamp;
            
            if isfield(obj.filteredData, 'altitude')
                alt = obj.filteredData.altitude;
            else
                alt = obj.data(1).attitude(:,4);
            end
            
            if nargin == 2
                figure(fig)
                hold on
            else
                fig = figure('name',sprintf('%s Altimeter Measurments',obj.name), 'numbertitle', 'off');
            end
            
            if numel(times) == numel(alt)
                plot(times, alt)
            else
                warning("Plot_Altitude --> Times and Alt are different sizes")
                plot(alt)
            end
            
            title(sprintf('%s: Atitude vs Time', obj.name), 'FontSize', 14)
            xlabel('Date-Time', 'FontSize', 14)
            ylabel('Atitude [m]', 'FontSize', 14)
            drawnow
            legend('FontSize', 12)
        end
        
        
        function fig = Plot_Alt_Geotiff(obj, geotiff)
            
            if isfield(obj.filteredData, 'altitude')
                alt = obj.filteredData.altitude;
            else
                alt = obj.data(1).attitude(:,4);
            end
            
            oo = alt == 0;
            
            lat_lon = obj.path.GT.lat_lon;
            
            fig = figure('name',sprintf("%s Altimeter", obj.name), 'numbertitle', 'off');
            geotiff.Show
            hold on
            
            plot(lat_lon(oo,2), lat_lon(oo,1), 'r.')        % Plot 0s
            plot(lat_lon(~oo,2), lat_lon(~oo,1), 'g.')      % Plot ~0s
            
            hold off
            
        end
        
        
        function fig = Show_Manifest(obj)
            
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
            % %             uit.Units =  'Normalized';
            %             uit.Position = [20, 20, 100, 100];
            
        end
        
        
        function fig = Plot_TimeStamps(obj)
            
            times = obj.data.timeStamp;
            
            fig = figure('name','Mission Timestamps','numbertitle','off');
            plot(times)
            
            drawnow
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



