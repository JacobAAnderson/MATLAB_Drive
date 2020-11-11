% Riptide_AcommsLog2Mat
% Jacob Anderson
% Robotic Decition Making Laboratory (RDML)
% Oregon State University
% Corallis OR
% August 20, 2020


% EcoMapperLog2Mat  ==========================================================================================================================================================
function Data = Riptide_AcommsLog2Mat(file, filters)

if nargin == 2
    rawData = Riptide_Read_Log(file, filters);                        % Read the log with filters
else
    rawData = Riptide_Read_Log(file);                                 % Read the log without filters
end


% Pars Acomms data ----------------------------------------------------------------------------------------
recive = rawData.ACOMMS_RECV_CSV;
sent   = rawData.ACOMMS_XMIT;
src    = rawData.src;
dest   = rawData.dest;
msg    = rawData.msg;


format = 'yyyy-MM-dd  HH:mm:ss.SSSSSS';
mission_cells = rawData.Mission;

for ii = numel(mission_cells): -1: 1
    
    % Pars sent message ---------------------------------------------------
    if isempty(sent{ii})
        sent_msg(ii)   = NaT('Format', format);
        sent_typ(ii,:) = ' ';
    else
        parts = strsplit(sent{ii}, ';');
        num   = uint64(str2num(parts{1}));
        sent_msg(ii)   = datetime(num, 'convertfrom', 'epochtime','TicksPerSecond',1e6, 'Format', format);
        sent_typ(ii,:) = parts{2};
    end
    
    
    
    % Pars recived message ------------------------------------------------
    if isempty(recive{ii})                  % No message recived
        recive_time(ii)  = NaT('Format', format);
        recive_msg(ii)   = NaT('Format', format);
        recive_typ(ii,:) = ' ';
        src_dest{ii}     = '\0';
    
    else                                    % Message Recived
        % Timestamp
        parts = strsplit(recive{ii}, '=');
        num   = str2double(parts{2});
        recive_time(ii) = datetime(num, 'convertfrom', 'epochtime', 'Format', format);

        % Source - Destination
        parts        = strsplit(src{ii}, '=');
        src_dest{ii} = parts{2};
        parts        = strsplit(dest{ii}, '=');
        src_dest{ii} = src_dest{ii} + "-" + parts{2};
 
        % Message
        parts = strsplit(msg{ii}, '=');
        parts = strsplit(parts{2}, ';');
        num   = uint64(str2num(parts{1}));
        recive_msg(ii)   = datetime(num, 'convertfrom', 'epochtime','TicksPerSecond',1e6, 'Format', format);
        recive_typ(ii,:) = parts{2};
    end
    
end


% Assign Values to Data Structure --------------------------------------------------------------------------

gpsLat = rawData.('GPS_LATITUDE');
gpsLon = rawData.('GPS_LONGITUDE');
depth  = rawData.DEPTH;


Data.header = {'  timeStamp: Date and time [MM.dd.yyyy HH:mm:ss.SS]';
               '    vehicle: Latitude [deg.dd], Longitude [deg.dd], depth from surface [m]';
               '     header: Explains the entries in this data structure';
               '        log: Name of the .log file that the data came from'};

Data.log       = rawData.log;
Data.timeErr   = rawData.timeErr;
Data.timeStamp = rawData.Sys_Time;
Data.gpsDate   = rawData.GPS_dateTime;
Data.mission   = rawData.Mission;
Data.modemID   = rawData.MODEM_ID;

Data.vehicle      = [gpsLat, gpsLon, depth];
Data.sent_msg     = sent_msg';
Data.sent_type    = sent_typ;
Data.src_dest     = src_dest';
Data.recived_msg  = recive_msg';
Data.recived_type = recive_typ;

Data.pAcomms_tof  = NaT(numel(rawData.Mission),4, 'Format', format);

Data.tof          = seconds( NaN(numel(rawData.Mission),1) );
end




