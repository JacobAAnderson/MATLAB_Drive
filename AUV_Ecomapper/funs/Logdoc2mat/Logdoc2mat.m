% Logdoc2mat
% Jacob Anderson
% Dec 10, 2018
% anderja7@oregonstate.edu

% Matlab function to import Starfish sidescan sonar binary files ".logdoc"
% Information on the Starfish 453 OEM hardware, which this script is customized for, can be found at:
%   https://www.blueprintsubsea.com/pages/product.php?PN=BP00755



function Data = Logdoc2mat(file)

[~, fileName, ext] = fileparts(file);

fprintf(" * Importing Side Scan Sonar File: %s\n", fileName)

%% Create Data Structure
Data = struct('timeStamp', datetime, 'position',  [],     'utmZone',   '',     ...
              'heading',   [],       'depth',     [],     'speed',     [],     ...
              'satellite', [],       'portSonar', struct, 'starSonar', struct, ... 
              'header',    [],       'settings',  '',     'log',       ''  );
          
Data.header = [ "timeStamp: dd/MM/yyyy hh:mm.SS.ssss";
                " position: Vehicle's position in utms [ Easting, Northing, Altetude ]";
                "  utmZone: Longitude and Latitude bands";
                "  heading: Vehicle's heading [degrees]";
                "    depth: Water depth [meters]";
                "    speed: Vehicle's speed [meters per second]";
                "portSonar: Sonar data structure";
                "starSonar: Sonar data structure";
                "      log: Name of the .log file that the data came from";
                "_______________________________________________________________________________________________________________";
                "Sonar data structure --> contrast: Sonar's 'Contrast' control setting  in decibels";
                "Sonar data structure --> offset:   Sonar's 'Offset'   control setting  in decibels ";
                "Sonar data structure --> gain:     Sonar's 'Gain'     control setting  in decibels";
                "Sonar data structure --> range:    Sonar's 'Range'    control setting  in Meters";
                "Sonar data structure --> signal:   Signal intensities of acoustic data in decibles ";
                "Sonar data structure --> EchoStrength: echo strength for visulization. Nomalized value, may contain NaN values";
                "Sonar data structure --> SamplePeriod: Sample period for the individual acoustic pings in seconds";
                "Sonar data structure --> Resolution:   Range Resolution for each sample in meters"];

Data.log = [fileName,ext];    % Add name of the log file to the data structure
Data.timeStamp = NaT;           % Default time value incase the file doesnt load

% Check that the file path is valid
if ~exist(file,'file')
    fprintf('\tThe File %s Does Not Exist\n', fileName);
    return
end


% Check that the files has data. i.e.. more than just a header
logfile = dir(fullfile(file));

if (logfile.bytes < 1000 )
    fprintf('\t--> File is empty\n')
    return
end

%% Read In File
fileID = fopen(file);
A = fread(fileID);
fclose(fileID);

A = uint8(A);


%%  Reader Header Block of Data File
% fileIdentifier = strtrim( convertCharsToStrings( char(A(1:32)) ));
% fileVersion    = typecast(A(33:36), 'uint32');
settingsLenght = typecast(A(37:40), 'uint32');
settingsStartAddress = 129;


%% Read Settings Block of Data File

VOS = 1475; % Speed of sound in water [m/s]

Data.settings = ParseSettings( A(settingsStartAddress : settingsStartAddress + settingsLenght));



%% Read Data Block of Data File

dataStartAddress = settingsStartAddress + settingsLenght; % Starting address of the first data message

% Prealocate Memeory
preSize = 10000;

bytes = uint8(1:200);

heading(preSize)      = 0;        
position(preSize,1:3) = [0, 0, 0];
utmZone(preSize)      = " ";
speed(preSize)        = 0;
sat(preSize)          = 0;
depth(preSize)        = 0;
dateTime(preSize)     = datetime;
portSonar(preSize)    = Scanline(bytes,VOS);
starSonar(preSize)    = Scanline(bytes,VOS);

index = 1;

% Cycle through the rest of the the data file and extract info
while dataStartAddress < length(A)
    
    % Make sure that the algoruthim is in sync with th elog file -------------------------------------------------------------------
    Sync = typecast(( A(dataStartAddress : dataStartAddress + 3) ), 'uint32');
    
    if ~ isequal(Sync, hex2dec('0FFFFFFFF'))
        disp('!!!!!!!!  Parsing in Not Syncronized   !!!!!!!!')
        return
    end
    
    % Pars payload data and metadata -----------------------------------------------------------------------------------------------
%   messageFormat = A( dataStartAddress + 4 );
    payloadType   = A( dataStartAddress + 5 );
    timeStamp     = Daytime( A( dataStartAddress + 6 : dataStartAddress + 13 ) );
%   sourceID      = A( dataStartAddress + 14 );
%   payloadFormat = A( dataStartAddress + 15 );
    payloadLength = typecast( A( dataStartAddress + 16 : dataStartAddress + 19 ), 'uint32');
    payload       = A( dataStartAddress+20 : dataStartAddress + 19 + payloadLength );
    checksum      = dec2bin(A( dataStartAddress : dataStartAddress + 19 + payloadLength ));
    
    % Do checksum -------------------------------------------------------------------------------------------------------------------
    [row,~] = size(checksum);
    check = xor(checksum(1,:), checksum(2,:));
    for ii = 3: row
        check = xor(check, checksum(ii,:));
    end
    
    if any(check)
        disp('!! ++++++  I no checksum --------- !!!')
    end
    
    % Interpret payload and assign to data structure -------------------------------------------------------------------------------------------------------
    
    switch dec2hex(payloadType)
        
        case '00' % Unknown Data type
            disp(' Unknown Data type')
            
            
        case '10' % Raw Data Message
            disp('Raw Data type')
            
            
        case '20' % Heading Message
            % disp('Heading Mesage')
            heading(index) = typecast(payload, 'double');                   % Heading in degrees as a decimal
            
            
        case '21' % Position Messages
            % disp('Position Mesage')
            [XX, YY, ZZ, zone] = PositionUTM(payload);                          % [Easting, Northing, Altetude, UTM zone] in meters
            position(index,1:3) = [XX, YY, ZZ];
            utmZone(index)  = zone;
            
            
        case '22' % Velocity Message
            % disp('Velocity Measuremnt')
            speed(index) = typecast(payload, 'double') * 0.27778;          % Velocity in meters per second
            
            
        case '23' % Satellite Message
            % disp('Satelite Message')
            sat(index) = Satellite(payload);                                   % Satellite fix and Dilution Of Position info
            
        case '24' % Depth Message
            % disp('Depth Message')
            depth(index) = typecast(payload, 'double');                      % Depth in meters
            
 
        case '26' % Sidescan Message
            % disp('Sidescan Message')
            dateTime(index) = timeStamp;
            
            portSonar(index) = Scanline( payload( 1 : end/2 ),    VOS ); % Sonar data structure
            starSonar(index) = Scanline( payload( end/2+1 : end), VOS );
            
            index = index +1;
            
        case '29'
            % I dont know what this is but it happens --> look into it later
            
        otherwise
            fprintf('\n  ! Unrecognized Payload Type: %d\n', payloadType )
           
            if payloadLength == 8 
                disp( typecast(payload, 'double') )
            end
    end
    
    dataStartAddress = dataStartAddress + 21 + payloadLength;                   % Calculate starting adress of the next data message
end


% Clean up un-ussed space -----------------------------------------------------
heading(index:preSize)    = [];        
position(index:preSize,:) = [];
utmZone(index:preSize)    = [];
speed(index:preSize)      = [];
sat(index:preSize)        = [];
depth(index:preSize)      = [];
dateTime(index:preSize)   = [];
portSonar(index:preSize)  = [];
starSonar(index:preSize)  = [];


% Clean up missing data -------------------------------------------------------
missingData = ismissing(utmZone);

if any(missingData)
    
    heading(missingData)    = [];       
    position(missingData,:) = [];
    utmZone(missingData)    = [];
    speed(missingData)      = [];
    sat(missingData)        = [];
    depth(missingData)      = [];
    dateTime(missingData)   = [];
    portSonar(missingData)  = [];
    starSonar(missingData)  = [];
    
end

% Assign Data to Data Struct --------------------------------------------------
Data.heading   = heading';
Data.position  = position;
Data.utmZone   = char(utmZone');
Data.speed     = speed';
Data.satellite = sat';
Data.depth     = depth';
Data.timeStamp = dateTime';
Data.portSonar = portSonar';
Data.starSonar = starSonar';

end



%% Functions ==================================================================================================================
function settings = ParseSettings(bytes)
bytes = bytes';
% disp('Sonar Settings')

settings = struct('Environment', '', 'Units', '', 'Log', '', 'Hardware', '');

string = char(bytes);

ind_uni = strfind(string,'\Units');
ind_env = strfind(string,'\Environment');
ind_log = strfind(string,'\Log');
ind_har = strfind(string,'\Hardware');
ind_fil = strfind(string,'\Filter');

ind_uni = [ind_uni, ind_env(1)];
ind_env = [ind_env, ind_log(1)];
ind_log = [ind_log, ind_har(1)];
ind_har = [ind_har, ind_fil(1)];
ind_fil = [ind_fil, numel(bytes)];

% Parse Units
for ii = 1: numel(ind_uni)-1
   str = strsplit( string(ind_uni(ii):ind_uni(ii+1)-1), char(0) );
   bit = bytes(ind_uni(ii):ind_uni(ii+1)-1);
   ind = find(bit == 0);
   settings.Units.(str{2}) = SettingValues(bit(ind(2)+1:end)); 
end

% Parse Environment
for ii = 1: numel(ind_env)-1
   str = strsplit( string(ind_env(ii):ind_env(ii+1)-1), char(0) );
   bit = bytes(ind_env(ii):ind_env(ii+1)-1);
   ind = find(bit == 0);
   settings.Environment.(str{2}) = SettingValues(bit(ind(2)+1:end));  
end

% Parse Log
for ii = 1: numel(ind_log)-1
   str = strsplit( string(ind_log(ii):ind_log(ii+1)-1), char(0) );
   bit = bytes(ind_log(ii):ind_log(ii+1)-1);
   ind = find(bit == 0);
   settings.Log.(str{2}) = SettingValues(bit(ind(2)+1:end)); 
end

% Parse Hardware
for ii = 1: numel(ind_har)-1
    bit = bytes(ind_har(ii):ind_har(ii+1)-1);
    str = strsplit( char(bit), char(0));
    str2 = strsplit( str{1}, {'\','[',']'});
    str3 = strsplit(str{2},{'[',']'});
    
    if size(str3,2) > 1, str3{1} = [str3{1},'_',str3{2}]; end
        
    ind = find(bit == 0);
    
    if size(str2,2) == 2
        continue        
    else
        settings.(str2{2}).(str2{3})( round(str2double(str2{4}))+1 ).(str3{1}) = SettingValues(bit(ind(2)+1:end));
    end
    
end

% Parse filter
% for ii = 1: numel(ind_fil)-1
%     bit = bytes(ind_fil(ii):ind_fil(ii+1)-1);
%     char(bit)
%     str = strsplit( char(bit), char(0))
%     
%     if size(str,2) == 1, continue, end
%     
%     str2 = strsplit( str{1}, {'\','[',']'})
%     str3 = strsplit(str{2},{'[',']'})
%     
%     if size(str3,2) > 1, str3{1} = [str3{1},'_',str3{2}]; end
%         
%     ind = find(bit == 0);
%     
%     settings.(str2{2}).(str3{1}) = SettingValues(bit(ind(2)+1:end));
%     
% end


end


function value = SettingValues(bytes)

type = typecast(bytes(1:4), 'uint32');

switch type
    case 0
        dips('LogDoc: Unsupported Data Type')
        
    case 1
        value = char(bytes(5:end));
        
    case 2
        value = bytes(5:end);
        
    case 4
        value = typecast(bytes(5:end), 'uint32');
        
    case 7
        value = bytes(5:end);
        
    case 11
        value = typecast(bytes(5:end), 'double');
        
    otherwise
        disp('LogDoc: Data format Unreckognized')
end
end


% conver Values -----------------------------------------------------------
function [valueOut, type] = convertValue(valueType,valueIn)

switch valueType
    case 'n/a'
        valueOut = NaN;
        type = '%f';
        
    case 'stringUtf16'
        valueOut = char(valueIn);
        type = '%s';
        
    case 'string-Array'
        valueOut = char(valueIn);
        type = '%s';
        
    otherwise
        valueOut = typecast(valueIn,valueType);
        type = '%d';
        
end
end


% Interpret timestamp data ------------------------------------------------
function timeStamp = Daytime(bytes)

% bytes(3:6) are [Month, Day, Hour, Minut] as a number representing their value.

year      = typecast( bytes(1:2), 'uint16');
seconds   = double(typecast( bytes(7:8), 'uint16'))/1000;

timeStamp = datetime([year, bytes(3:6)', seconds], 'Format', 'dd/MM/yyyy hh:mm.SS.ssss');
end


function [x,y,z, zone] = PositionUTM(bytes)                 % Interpret loaction data

zone = strcat(string(bytes(2)), {' '}, char(bytes(1)));     %  Longitude and Latitude bands (respectivly)

y = typecast(bytes( 3:10), 'double');   % Northing
x = typecast(bytes(11:18), 'double');   % Easting
z = typecast(bytes(19:26), 'double');   % Altitude

end


function  [fix, fix3D, PDOP, HDOP, VDOP] = Satellite(bytes) % Interpret satalite data

fix   = logic(bytes(1));                        % Indicates valid GPS fix
fix3D = logic(bytes(2));                        % Indicates 3D GPS fix
PDOP  = typecast(bytes( 3: 6), 'single');       % Positional Dilution Of Position of the fix
HDOP  = typecast(bytes( 7:10), 'single');       % Horizontal Dilution Of Position of the fix
VDOP  = typecast(bytes(11:14), 'single');       % Verticle   Dilution Of Position of the fix

end


function SS = Scanline(bytes,VOS)                           % Interpret Sonar data
% Pars Sonar data --------------------------------------------------------------------------------------------------------
SS.contrast = typecast(bytes( 1: 4), 'single');     % Sonar's "Contrast" control setting  in decibels
SS.offset   = typecast(bytes( 5: 8), 'single');     % Sonar's "Offset"   control setting  in decibels
SS.gain     = typecast(bytes( 9:12), 'single');     % Sonar's "Gain"     control setting  in decibels
SS.range    = typecast(bytes(13:16), 'single');     % Sonar's "Range"    control setting  in Meters
SS.signal   = single(bytes(17:end))./2.0;           % Signal intensities of acoustic data in decibles
                                                    %   * This values where origonaly multiplied by 2 to maximixe storage space.


% Interpret Sonar data ----------------------------------------------------------------------------------------------------
echoStrength = (SS.signal - SS.offset + SS.gain)./SS.contrast;  % Calculate echo strength for visulization --> Normalized
echoStrength(echoStrength > 1) = NaN;                           % Values must be between 0 and 1
echoStrength(echoStrength < 0) = NaN;

SS.EchoStrength = echoStrength;

samples         = length(bytes) - 16;           % Number of samples
duration        = SS.range*2/VOS;               % Duration of the scan
SS.SamplePeriod = duration/samples;             % Sample period for the individual acoustic pings in seconds
SS.Resolution   = SS.range/samples;             % Range Resolution for each sample in meters

end



