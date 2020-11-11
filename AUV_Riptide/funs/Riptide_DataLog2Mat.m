% Riptide_DataLog2Mat
% Jacob Anderson
% Robotic Decition Making Laboratory (RDML)
% Oregon State University
% Corallis OR
% August 20, 2020


% EcoMapperLog2Mat  ==========================================================================================================================================================
function Data = Riptide_DataLog2Mat(data_file, offset_file, filters)

if nargin == 3
    [rawData, sz] = Riptide_Read_Log(data_file, filters);                        % Read the log with filters
else
    [rawData, sz] = Riptide_Read_Log(data_file);                                 % Read the log without filters
end

% ---- Correct Bathymetry Data for the orientation of the vehicle ---------
date     = rawData.('GPS_dateTime');
gpsLat   = rawData.('GPS_LATITUDE');
gpsLon   = rawData.('GPS_LONGITUDE');
depth    = rawData.('DEPTH');
altitude = rawData.('ALT_ALTITUDE');
roll     = rawData.('NAV_ROLL');
pitch    = rawData.('NAV_PITCH');
heading  = rawData.('NAV_HEADING');


% Apply water height off set
if exist(offset_file, 'file')
    altitude = ApplyOffset(offset_file, altitude, date);    
end


[x,y,utmzone] = deg2utm(gpsLat,gpsLon);                                     % Vehilces position in utm --> [x,y,utmzone] = deg2utm(Lat,Lon)

riptide = num2cell([x,y,depth]',1);                                         % Vehicles position in 3 space
alt_pitch = pitch + 21.42;                                                  % Adjust for altimiter mounting angle

p_rt  = arrayfun( @(z) [0;0;z],                                                         altitude,           'UniformOutput',false); % Depth reading as a vector under the vehicle
Rx    = arrayfun( @(roll) [ 1 0 0; 0 cos(roll) -sin(roll); 0 sin(roll) cos(roll)],      -roll,              'UniformOutput',false); % X rotation matrix with role angel
Ry    = arrayfun( @(pitch)[ cos(pitch) 0 sin(pitch); 0 1 0; -sin(pitch) 0 cos(pitch)],  -alt_pitch,         'UniformOutput',false); % Y rotation matrix with pitch angle
Rz    = arrayfun( @(yaw)  [ cosd(yaw) -sind(yaw) 0; sind(yaw) cosd(yaw) 0; 0 0 1],      -heading,           'UniformOutput',false); % Z rotation matrix with compass heading
R     =  cellfun( @(a, b, c) a*b*c,                                                     Rz, Ry, Rx,         'UniformOutput',false); % Compounded rotation matrix
P_w   =  cellfun( @(r,eco,p) [ r, eco; 0 0 0 1] * [p;1],                                R, riptide', p_rt,  'UniformOutput',false); % Coordinate transormation

P_w = cell2mat(P_w);
P_w = reshape(P_w,[4,sz])';

[Lat_ ,Lon_] = utm2deg(P_w(:,1),P_w(:,2),utmzone);                          % Convert the bottom point's coordinates into lat lon



% Assign Values to Data Structure --------------------------------------------------------------------------
gpsLat = rawData.('GPS_LATITUDE');
gpsLon = rawData.('GPS_LONGITUDE');
depth  = rawData.DEPTH;
speed  = rawData.('RT_THRUST_SPEED');


Data.header = {'  timeStamp: Date and time [MM.dd.yyyy HH:mm:ss.SS]';
               '    vehicle: Latitude [deg.dd], Longitude [deg.dd], depth from surface [m]';
               '      bathy: Latitude [deg.dd], Longitude [deg.dd], depth [m]';
               '   attitude: Roll [deg], Pitch [deg], Yaw [deg], Altitude [m], Speed [m/s]'
               '     header: Explains the entries in this data structure';
               '        log: Name of the .log file that the data came from';
               'Known Stats: Altimeter: 0 mean, 0.059687[m] std' };


Data.log       = rawData.log;
Data.timeErr   = rawData.timeErr;
Data.timeStamp = rawData.Sys_Time;
Data.gpsDate   = rawData.GPS_dateTime;
Data.mission   = rawData.Mission;

Data.vehicle    = [gpsLat, gpsLon, depth];
Data.bathymetry = [Lat_, Lon_, P_w(:,3)];
Data.attitude   = [roll, pitch, heading, altitude, speed];
Data.tide       = zeros(sz,1);
end




function alt = ApplyOffset(file, alt, date)

fileID = fopen(file,'r');
[~] = fgetl(fileID); 
offsets = textscan(fileID, '%{yyyy-MM-dd}D%f', 'Delimiter', ',', 'ReturnOnError', false); % read the rest of the file
fclose(fileID);

for valid = find(~isnat(date))                                                 % Find a Valid time reading
    
    date = date(valid(1));
    date = datetime(year(date), month(date), day(date));
    
    in = ismember(offsets{1}, date);
    
    if ~any(in)
        warning('Mission Date does not correspond to any dates in the off set log')
        return
    end
    
    dif = offsets{2};
    dif = dif(in);
    
    oo = alt ~= 0;
    
    alt(oo,:) = alt(oo,:) - dif;
    
    fprintf("\tApplying Altitude Offest of: %f\n", dif)
    
    return
end

end







