

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

