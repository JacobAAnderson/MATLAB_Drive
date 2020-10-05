% Read Riptide file
% Jacob Anderson
% Sept 10, 2020


function [rawData, sz] = Riptide_Read_Log(file, filters)

[~, fileName, ext] = fileparts(file);

fprintf(' * Importing Riptide Data Log: %s\n', [fileName,ext]);

            
rawData.log = [fileName,ext];    % Add name of the log file to the data structure

rawData.timeStamp = NaT;         % Default timestamp value incase the file is unreadable

% Check that the file path is valid
if ~exist(file,'file')
    fprintf('\tThe File %s Does Not Exist\n', fileName);
    return
end


% Check that the files has data. i.e.. more than just a header
logfile = dir(file);
if (logfile.bytes < 1000 )
    fprintf('\t--> File is empty\n')
    return
end


% Read txt file
try 
    fileID = fopen(file,'r');
    Names = strsplit(fgetl(fileID), ',');                                   % Get header to use as variable names
    Names = strtrim(Names);                                                 % Get rid of white spaces
    formatSpec = Riptide_LogFormat(Names);                                  % Determin data format based on variable
    dataArray = textscan(fileID, formatSpec, 'Delimiter', ',', 'ReturnOnError', false); % read the rest of the file
    fclose(fileID);
    
catch exception                 % File Cannot be read --> Skip it and return empty struct
    warning off backtrace
    warning('Error Reading log file \n\t\t%s \n\t\t%s \n\t\tError in: %s, Line: %d \n\n\t\tSuspect bad log file, skipping it \n\n',exception.identifier,exception.message,exception.stack(1).name,exception.stack(1).line)
    warning on backtrace
    return
    
end


% ---- Check for errors in the Acomms log ---------------------------------
% --- Fill in empty lines ---
L = length(dataArray{1});
for ii = 1: numel(dataArray)
    if L > length(dataArray{ii})
        warning('Filling in empty lines in row %d of %d rows',ii, numel(dataArray))
        a = dataArray{ii}; 
        a{L} = ""; 
        dataArray{ii} = a;
    end
end


% ---- Data filters Data by GPS fix ----------------
in = true(L,1);
if nargin == 2
    for filter = filters'
        
        tf = strcmp(Names, filter{1});
        
        if any(tf)
            in(dataArray{tf} == filter{2} )   = false;
        end
        
    end
end

if nargout ==2, sz = sum(in); end

if ~any(in), warning('This Data Log has no usable data'), end               % All the data has been filtered out


% ---- Extrace the data from their cell arrays ----------------------------
for ii = 1:numel(Names)
    data = dataArray{ii};
    rawData.(Names{ii}) = data(in);                                         % Assign data to struc using variable names as fields
end

rawData.DEPTH = zeros(sum(in),1);                                           % Depth data is currently not availalbe



% -----Get the mission names seperated from the paths ---------------------
mission = rawData.('Mission');

tf = cellfun(@isempty, mission);

mission(tf) = {'idel'};

% Determin if speed data is available
if isfield(rawData, 'RT_THRUST_SPEED')                                      % Speed data is available for Data logs but not all Acomms
    speed = rawData.('RT_THRUST_SPEED');
    
    tf2 = speed == 0;
    mission(tf2) = {'idel'};                                            % If Speed is 0, then the vehicle is idel
    
    tf = tf | tf2;
end


for ii = numel(mission): -1: 1
    
    if ~tf(ii)                                    
        parts = strsplit(mission{ii}, '/');                           % Get Mission name out of text string
        mission{ii} = parts{5};
    end
end

rawData.('Mission') = mission;



% Compensate for bad time stamp -----------------------------------------------------------------------------
dates = rawData.('GPS_dateTime');
if isfield(rawData, 'Sys_Time_')                                            % Date and time of the measurment
    
    times = rawData.('Sys_Time_');
    
    timestep = [0; times(2:end) - times(1: end-1)];                         % Calculate step between time
    
    for f = find( timestep < 0 )                                            % Find role overs where time step is negative
        fprintf("\tCorrecting For Role-over\n")
        times(f:end) = times(f:end) + hours(1);                             % Add an hour onto the subsiquent time entries to compensate for the role over              
    end
    
    timestep = [0; times(2:end) - times(1: end-1)];                         % Re-calc the time steps
    
    timeStamp = datetime(dates, 'InputFormat', 'yy/MM/dd HH:mm:ss', 'Format', 'yyyy-MM-dd  HH:mm:ss.SSSSSS');
    
    timeStamp  = timeStamp(1) + cumsum(timestep);                           % Use the frist GPS date and time to start the timestamps, and add time steps to that gps point
     
else, timeStamp = rawData.('Sys_Time');
end


% Try to sycronize timeStamp to GPS Time ----------------------------------
%  * For logs from berfor the system clocls were syncronized

if timeStamp(1) < datetime(2020,9,11)
    
     fprintf("\tTrying to sync system time to GPS time\n")
    
    diff = duration(0, 10, 0, 0);
    evalTime = timeStamp;
    
    timestep1 = [0; timeStamp(2:end) - timeStamp(1: end-1)];
    
    count = 1;
    total_Correction = duration(0,0,0,000);
    while abs(diff) >= duration(0,0,0,500) && count < 100
        
        err = evalTime - dates;
        
        diff = mean(err);
        
        evalTime = evalTime - diff;
        
        count = count + 1;
        
        total_Correction = total_Correction + diff;
    end
    
    timeStamp = evalTime;
    
    timestep2 = [0; timeStamp(2:end) - timeStamp(1: end-1)];
    
    cumErr = sum(abs(timestep1 - timestep2));
    cumErr.Format = 's';
    fprintf("\tCumulative Error:")
    disp(cumErr)
    
    if cumErr >= seconds(0.5)
        warning('Cumulative Error is clock sycronization is %D', cumErr)
    end

else
   total_Correction = 0; 
end

rawData.('Sys_Time') = timeStamp;
rawData.timeErr      = total_Correction;

end






function frmt_str = Riptide_LogFormat(header)

frmt_str = "";

for ii = 1:numel(header)
    switch header{ii}
        
        case 'Sys_Time_', f = '%{mm:ss.S}D'; 
        case 'Sys_Time',  f = '%{yyyy-MM-dd HH:mm:ss.SSSSSS}D';
                
        case 'Mission', f = '%s';
        
        case 'GPS_fix',       f = '%d';
        case 'GPS_dateTime',  f = '%{yy/MM/dd HH:mm:ss}D';
        case 'GPS_LATITUDE',  f = '%f';
        case 'GPS_LONGITUDE', f = '%f';
        
        case 'NAV_LAT',     f = '%f';
        case 'NAV_LONG',    f = '%f';
        case 'NAV_HEADING', f = '%f';
        case 'NAV_ROLL',    f = '%f';
        case 'NAV_PITCH',   f = '%f';
        case 'NAV_YAW',     f = '%f';
        
        case 'PWR_NOSE_VOLTAGE',    f = '%f';
        case 'PWR_TAIL_VOLTAGE',    f = '%f';
        case 'PWR_PAYLOAD_VOLTAGE', f = '%f';
        case 'PWR_NOSE_FAULT',      f = '%d';
        case 'PWR_TAIL_FAULT',      f = '%d';
        case 'PWR_PAYLOAD_FAULT',   f = '%d';
            
        case 'RT_THRUST_SPEED', f = '%f';
            
        case 'ALT_TRIGGER',     f = '%d';
        case 'ALT_PING_RATE',   f = '%d';
        case 'ALT_SOUND_SPEED', f = '%d';
        case 'ALT_ALTITUDE',    f = '%f';
        
        case 'MODEM_ID',	    f = '%d';
        case 'ACOMMS_XMIT',	    f = '%s';
        case 'ACOMMS_RECV_CSV', f = '%s';
        case 'src',             f = '%s';
        case 'dest',            f = '%s';
        case 'msg',             f = '%s';
            
        otherwise, f = '%s';
    end
    
    frmt_str = frmt_str +  f;
     
end

frmt_str = frmt_str + '%[^\n\r]';

end

