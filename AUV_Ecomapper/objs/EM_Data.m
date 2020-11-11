% Ecomapper Data Class
% Jacob Anderson
% 8/28/2019

classdef EM_Data
    
    properties
        RawData
        FilteredData
        Path
    end
    
    methods
        
        % Instanciate Class by loading Raw Data
        function obj = EM_Data(file), warning on backtrace
            
            if isfile(file)
                
                obj.RawData = EcoMapperLog2Mat(file,'dvl filter');
                
            elseif isfolder(file)
                
                fileType = 'log';
                fun = @(file_) EcoMapperLog2Mat(file_,'dvl filter');
                
                obj.RawData = BatchDir(file,fileType,fun);
                
            else
                disp("EM_Data --> Invalid File Path")
                return
            end
      
        end
        
        
        % Low Pass Filter Altimiter Data
        function obj = FilterAltimiter(obj, beta)
            
            for jj = 1: size(obj.RawData,1)
                
                data_size = size(obj.RawData(jj).attitude,1);
                
                obj.FilteredData(jj).Altimeter = zeros(data_size,1);
                
                evalWindow = SlidingWindow(10, data_size);
                adjustWindow = SlidingWindow(4, data_size);
                
                for ii = 1: data_size
                    
                    eval_indx = evalWindow.D1(ii);
                    eval_indx(eval_indx == ii) = [];
                    
                    data = obj.RawData(jj).attitude(eval_indx,4);
                    
                    datum = obj.RawData(jj).attitude(ii,4);
                    sigma = std(data);
                    mu    = mean(data);
                    
                    if abs(datum - mu) > beta * sigma
                        adjust_indx = adjustWindow.D1(ii);
                        adjust_indx(adjust_indx==ii) = [];
                        
                        obj.FilteredData(jj).Altimeter(ii) = mean(obj.RawData(jj).attitude(adjust_indx,4));
                    else
                        obj.FilteredData(jj).Altimeter(ii) = datum;
                    end
                    
                end
            end
        end
        
        
        % Low Pass Filter for Headings
        function obj = FilterHeading(obj, beta)
            
            % figure('name','Before')
            % plot( EM_Data.attitude(:,3))
            
            data_size = size(obj.RawData.attitude,1);
            
            evalWindow   = SlidingWindow(10, data_size);
            adjustWindow = SlidingWindow( 4, data_size);
            
            for ii = 1: data_size
                
                eval_indx = evalWindow.D1(ii);
                eval_indx(eval_indx == ii) = [];
                
                data = obj.RawData.attitude(eval_indx,3);
                
                delta = abs(data(2:end) - data(1:end-1));
                
                datum = obj.RawData.attitude(ii,3);
                
                
                adjust_indx = adjustWindow.D1(ii);
                adjust_indx(adjust_indx==ii) = [];
                
                adjust_data = obj.RawData.attitude(adjust_indx,3);
                
                if max(delta) > 180
                    in = data < 180;
                    data(in) = data(in) + 360;
                    datum = datum + 360;
                    
                    in = adjust_data < 180;
                    adjust_data(in) = adjust_data(in) + 360;
                end
                
                sigma = std(data);
                mu    = mean(data);
                
                if abs(datum - mu) > beta * sigma
                    
                    new_theta = mean(adjust_data);
                    
                    if new_theta > 360, new_theta = new_theta-360; end
                    
                    obj.FilteredData.Heading(ii) = new_theta;
                else
                    obj.FilteredData.Heading(ii) = obj.RawData.attitude(ii,3);
                end
            end
            
            % figure('name','After')
            % plot( EM_Data.attitude(:,3))
        end
        
        
        % Corect Depth Readings for Tides
        function obj = TideCorrection(obj, stationID, plotResults)
            
            if nargin < 3
                plotResults = false;
            end
            
            
            % Get Tide Data -----------------------------------------------------------
            times = cat(1, obj.RawData.timeStamp);
            
            start = datetime(times(1));         % Beging of time duration
            stop  = datetime(times(end));       % End of time duration
            
            tide = NOAA_tides( stationID, start, stop, plotResults);
            
            
            % Coralate with endata ----------------------------------------------------
            for jj = 1:size(obj.RawData,1)
                
                for ii = 1: numel(obj.RawData(jj).timeStamp)
                    
                    timeDiff = abs( tide.DateTime - obj.RawData(jj).timeStamp(ii) );
                    
                    indx = timeDiff == min(timeDiff(:));
                    
                    obj.RawData(jj).tide(ii) = tide.Height(indx);
                    
                end
                
                if plotResults
                    figure('Name','Get Tide Data - Results', 'NumberTitle','off')
                    plot(obj.RawData(jj).timeStamp, emData(jj).tide)
                    xlabel('Date-Time')
                    ylabel('Tide Height [m]')
                end
            end
        end
        
        
        % Get Vehicle Path in UTM an dlat-lon
        function obj = GetVehiclePaths(obj)
            
            for ii = 1: size(obj.RawData,2)
                
                speed     = obj.RawData(ii).attitude(:,5);
                heading   = obj.RawData(ii).attitude(:,3);
                timestamp = obj.RawData(ii).timeStamp;
                
                timeStep = seconds(timestamp(2:end) - timestamp(1: end-1));
                
                deltaX_dr = speed(1:end-1) .* timeStep .* sind( heading(1:end-1) );
                deltaY_dr = speed(1:end-1) .* timeStep .* cosd( heading(1:end-1) );
                
                x = [0; cumsum(deltaX_dr)];
                y = [0; cumsum(deltaY_dr)];
                
                [xx, yy, utmzone] = deg2utm( obj.RawData(ii).vehicle(:,1), obj.RawData(ii).vehicle(:,2) );
                
                [dr_lat, dr_lon] = utm2deg(x + xx(1), y+yy(1), utmzone);
                
                obj.Path(ii).GT.lat_lon = obj.RawData(ii).vehicle(:,1:2);
                obj.Path(ii).GT.utm     = [xx,yy];
                obj.Path(ii).GT.utmzone = utmzone;
                
                obj.Path(ii).DR.utm(:,1) = x + xx(1);
                obj.Path(ii).DR.utm(:,2) = y + yy(1);
                obj.Path(ii).DR.lat_lon  = [dr_lat, dr_lon];
                
            end
        end
        
        
        % Get DR Path Structure
        function obj = Altimeter_Path(obj)
            
            for ii = 1: size(obj.RawData,2)
                
                [xx, yy, utmzone] = deg2utm(obj.RawData.bathymetry(:,1), obj.RawData.bathymetry(:,2));
                
                path_diff = obj.Path(ii).DR.utm - obj.Path(ii).GT.utm;
                
                obj.Path(ii).alt.lat_lon = obj.RawData.bathymetry(:,1:2);
                obj.Path(ii).alt.utm     = [xx, yy];
                obj.Path(ii).alt.utmzone = utmzone;
                obj.Path(ii).alt.DR      = obj.Path(ii).alt.utm + path_diff;
                
            end
        end
        
        
        % Plotting
        function PlotPaths(obj, Geotiff)
            
            ii = 1;
            
            figure('name', 'Ecomapper Pathes', 'numbertitle', 'off')
            
            if nargin > 1, geoshow(Geotiff.Image, Geotiff.Data), end
            
            hold on
            
            p1 = plot(obj.Path(ii).GT.lat_lon(:,2), obj.Path(ii).GT.lat_lon(:,1), 'g');     % Plot Pathe
            p2 = plot(obj.Path(ii).DR.lat_lon(:,2), obj.Path(ii).DR.lat_lon(:,1), 'r');
            
            plot(obj.Path(ii).GT.lat_lon(1,2), obj.Path(ii).GT.lat_lon(1,1), 'dm');         % Plot Begining
            
            plot(obj.Path(ii).DR.lat_lon(end,2), obj.Path(ii).DR.lat_lon(end,1), '*r');     % Plot Ends
            plot(obj.Path(ii).DR.lat_lon(end,2), obj.Path(ii).DR.lat_lon(end,1), 'or');
            
            plot(obj.Path(ii).GT.lat_lon(end,2), obj.Path(ii).GT.lat_lon(end,1), '*g');
            plot(obj.Path(ii).GT.lat_lon(end,2), obj.Path(ii).GT.lat_lon(end,1), 'og');
            
            hold off
            
            title('Ecomapper Paths')
            xlabel('Easitng')
            ylabel('Northing')
            legend([p1,p2],{'GT','DR'})
            
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
