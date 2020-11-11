% Function to Import data from a MOOS alog


function Data = MOOS_Alog2mat( fileName, filePath )

% Read file
file = fullfile(filePath, fileName);
fileID  = fopen(file,'r');
tline = 'go';


% Create Data Structure
Data = struct('gps',{},'nav',{},'alt',{},'Header','','Logs','');
Data(1).Header = 'gps: Date, Lon, Lat, X, Y, Heading, Speed, Fix, Quality';
Data(2).Header = 'nav: Pitch, Roll, Heading';
Data(3).Header = 'alt: Range';

Data(1).Logs = fileName;

navIndex  = 0;
gpsIndex  = 0;
altIndex  = 0;

thisTime  = [];
navTime(20) = 0;
altTime(20) = 0;
reCalNav = false;
reCalAlt = false;

gpsDate   = datetime();
gpsYear   = [];
gpsMonth  = [];
gpsDay    = [];
gpsHour   = [];
gpsMinute = [];

%for ii = 1:63
while ischar(tline)
    
    tline = fgetl(fileID);                                                  % Get the next line of the file
    
    if ~ischar(tline) || ~isequal(regexp(tline, '^\d*'), 1)                 % Make shure this line of data is something that we want
        continue
    end
    
    % '[A-Z]+\_[A-Z]+'
    % '\<GPS+\_[A-Z]+'
    
    
    if isequal(regexp(tline, '\<GPS+\_[A-Z]+', 'once'), 17)                 % Do stuff with GPS data
        
        data = strsplit( tline, ' ');
        time = str2double(data{1});
        
        switch data{2}
            
            case 'GPS_YEAR'
                gpsYear = 2000 + str2double(data{4});
                
            case 'GPS_MONTH'
                gpsMonth = str2double(data{4});
                
            case 'GPS_DAY'
                gpsDay = str2double(data{4});
                
            case 'GPS_HOUR'
                gpsHour = str2double(data{4});
                
            case 'GPS_MINUTE'
                gpsMinute = str2double(data{4});
                
            case 'GPS_SECOND'
                thisTime = time;                                            % Keep Track of time and clear the existing measurments
                gpsSecond = str2double(data{4}) + mod(abs(thisTime),1);
                gpsDate = datetime([gpsYear, gpsMonth, gpsDay, gpsHour, gpsMinute, gpsSecond], 'Format','MM/dd/yy HH:mm:SS.sss');
                
                Data(gpsIndex).gps{1} = gpsDate;
                
            case 'GPS_LONGITUDE'
                Data(gpsIndex).gps{2} = str2double(data{4});
                
            case 'GPS_LATITUDE'
                Data(gpsIndex).gps{3} = str2double(data{4});
                
            case 'GPS_X'
                Data(gpsIndex).gps{4} = str2double(data{4});
                
            case 'GPS_Y'
                Data(gpsIndex).gps{5} = str2double(data{4});
                
            case 'GPS_HEADING'
                Data(gpsIndex).gps{6} = str2double(data{4});
                
            case 'GPS_SPEED'
                Data(gpsIndex).gps{7} = str2double(data{4});
                
            case 'GPS_FIX'
                Data(gpsIndex).gps{8} = str2double(data{4});
                
            case 'GPS_QUALITY'                                              % This is first gps data measuremnt in the string of incoming gps data
                gpsIndex = gpsIndex + 1;
                
                Data(gpsIndex).gps{9} = str2double(data{4});
                
                gpsYear   = [];
                gpsMonth  = [];
                gpsDay    = [];
                gpsHour   = [];
                gpsMinute = [];
                
            otherwise
                continue
        end
        
        
    elseif isequal(regexp(tline, '\<NAV+\_[A-Z]+', 'once'), 17)             % Do stuff with Navigation data
        data = strsplit( tline, ' ');
        time = str2double(data{1});
        
        switch data{2}
            case 'NAV_PITCH'
                navIndex = navIndex + 1;
                
                if isempty(thisTime)
                    navTime(navIndex) = time;
                    
                elseif ~isempty(thisTime) && ~reCalNav
                    % disp('Calculate exisintg date-times')
                    reCalNav = true;
                    navTime(navIndex) = time;
                    
                    for ii = 1: navIndex
                        Data(ii).nav{1} = gpsDate + seconds( navTime(ii) - thisTime );
                    end
                    
                else
                    Data(navIndex).nav{1} = gpsDate + seconds( time - thisTime );
                    
                end
                
                Data(navIndex).nav{2} = str2double(data{4});
                
            case 'NAV_ROLL'
                Data(navIndex).nav{3} = str2double(data{4});
                
            case 'NAV_HEADING'
                Data(navIndex).nav{4} = str2double(data{4});
                
            otherwise
                continue
        end
        
        
    elseif isequal(regexp(tline, '\<ALT+\_[A-Z]+', 'once'), 17)             % Do stuff with Altimiter data
        data = strsplit( tline, ' ');
        time = str2double(data{1});
        
        switch data{2}
            case 'ALT_RANGE'
                altIndex = altIndex + 1;
                
                if isempty(thisTime)
                    altTime(altIndex) = time;
                    
                elseif ~isempty(thisTime) && ~reCalAlt
                    % disp('Calculate exisintg date-times')
                    reCalAlt = true;
                    altTime(altIndex) = time;
                    
                    for ii = 1: altIndex
                        Data(ii).alt{1} = gpsDate + seconds( altTime(ii) - thisTime );
                    end
                    
                else
                    Data(altIndex).alt{1} = gpsDate + seconds( time - thisTime );
                    
                end
                
                Data(altIndex).alt{2} = str2double(data{4});
                
            otherwise
                continue
        end
        
        
    else
        continue
        
    end
    
    
    
    
    
end

fclose(fileID);

disp('Finised Importing Alog')


end


