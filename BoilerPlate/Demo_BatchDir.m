%% Batch process a directory
close all
clear all
clc


%% Choose Data Type --> Uncomment the block that you want

%___ Load Ecomapper Mission Logs ________________________
fileType = 'log';
fun = @(file) EcoMapperLog2Mat(file,'dvl filter');

% %___ Load Ecomapper Side Scan Sonar Files ______________ 
% fileType = 'logdoc';
% fun = @(file) Logdoc2mat(file);



%% Load Data
ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(true);                                                     % Indicate whether new data should be selescted by the user

filePath = ui.GetDir('Data',fileType);                                      % Get folder GUI

if isempty(filePath)                                                        % Check if the Input box was cnaceled
    disp('Get Directory Input box Cancelled')
    return                                                                  % End script if there isn't an input for it to use
end

Data = BatchDir(filePath,fileType,fun);


if isempty(Data)                                                            % Check that data is present
    disp("Empty Data Structure")
    disp("Doulbe check the directroy that you chose")
    return                                                                  % End script if there isn't any data to process
end

disp('Data Structure:')
disp(Data(1).header)



%% Do Stuff With the Data
switch fileType
    
% Data from Ecomapper Mission Files ---------------------------------------------------------------
    case 'log'     
        
        % Get Data out of data structure
        vehicle         = cat(1, Data.vehicle);
        bathymetry      = cat(1, Data.bathymetry);
        waterParameters = cat(1, Data.wqData);
        
        lon   = vehicle(:,2);
        lat   = vehicle(:,1);
        depth = vehicle(:,3);
        
        [geoImage, geoData] = ui.getGeoTiff( mean(lon), mean(lat) );        % Gui to get Geotiff
        
        
        % Show GPS track ----------------------------------------------------------
        figure('Name','GPS track','NumberTitle','off')
        hold on
        if ~isempty(geoImage)                                                        % Check if the Data Input box was cnaceled
            geoshow(geoImage, geoData)
        end
        plot(lon,lat,'g')
        hold off
        
        
        % Plot some data ---------------------------------------------------------------------------
        temp  = waterParameters(:,9);                                       % Water Temperature [deg C]
        
        figure('Name','Water Temperature','NumberTitle','off')
        scatter3(lon,lat,-depth,'.','CData',temp)
        colorbar
        view(30,45)



% Data from Ecomapper Side Scan Sonar Files ------------------------------------------------------    
    case 'logdoc'
       
        % Get GPS Data for the option of downloading geotiff from the web.
        utm        = cat(1,Data(1).position);                               % Vehicle location in UTMs [ Easting, Northing, Altitude]
        utmzone    = cat(1,Data(1).utmZone);                                % UTM-zone: Longitude and Latitude bands
        [lat, lon] = utm2deg(utm(:,1), utm(:,2), utmzone);                  % Conver UTM to lat lon
        [geoImage, geoData] = ui.getGeoTiff( mean(lon), mean(lat) );        % Gui to get Geotiff

        
        % Cycle through the various scan lines
        for ii = 1: length(Data)
            
            % Get Data out of data structure
            utm       = cat(1,Data(ii).position);                           % Vehicle location in UTMs [ Easting, Northing, Altitude]
            utmzone   = cat(1,Data(ii).utmZone);                            % UTM-zone: Longitude and Latitude bands
            heading   = cat(1,Data(ii).heading);                            % Vehicle heading [deg]
            speed     = cat(1,Data(ii).speed);                              % Vehicle speed   [m/s]
            depth     = cat(1,Data(ii).depth);                              % Water depth     [m]
            portSonar = cat(1,Data(ii).portSonar);                          % Port Sonar data structure
            starSonar = cat(1,Data(ii).starSonar);                          % Starboard Sonar data structure
            
            
            % Display the vehicle's location when the data was collected -------------------------------------------------------------
            [lat, lon] = utm2deg(utm(:,1), utm(:,2), utmzone);              % Conver UTM to lat lon
                       
            if isempty(geoImage)                                            % Check if the Data Input box was cnaceled
                disp('Input box Cancelled')                                 % If no geotiff was selected just plot the location data
                figure('Name','GPS data from Sonar file','NumberTitle', 'off')
                hold on
                scatter(lon,lat,'.')
                p1 = plot(lon(1),lat(1),'*','Color','r');                   % Indicate the starting position
                legend(p1,'Starting Location')
                hold off
                
            else
                figure('Name','GPS data from Sonar file','NumberTitle', 'off') % Plot location data over geotiff
                geoshow(geoImage, geoData)
                hold on
                scatter(lon,lat,'.')
                p1 = plot(lon(1),lat(1),'*','Color','r');                   % Indicate the starting position
                legend(p1,'Starting Location')
                hold off
            end
            
            % Display Sonar image data -----------------------------------------------------------------------------------------------
            portEcho = cat(2,portSonar.EchoStrength);                       % Concatinate the image date from the sonar data structures
            starEcho = cat(2,starSonar.EchoStrength);
            
            image = [flipud(portEcho); starEcho];                           % Put the port and starboard images together
            
            figure('Name', 'Sonar Image','NumberTitle','off')               % Show the image
            hold on
            imshow( image )
            l1 = line([0,0],[0,size(image,1)],'color','r');
            legend(l1,'Beginning Of Sonar Track', 'Location','northwest')
            hold off
            
        end
        
end
