% Logdoc2mat
% Jacob Anderson
% Dec 10, 2018
% anderja7@oregonstate.edu

% Matlab function to import Starfish sidescan sonar binary files ".logdoc"

%clear all
close all
clc

UserInputs = SavedUserInputs(mfilename);                                    % Instantiate the Saved User Inputs Class for working with file paths

%% Get Sonar data
% Get file and folder GUI ----------------------------------------------------------------------------------------------
file = UserInputs.getFile('SonarData','sds');            % GUI to get the name and file path of a file
fprintf('\nFile Path and name: %s \n\n', file);

if isempty(file)                                                            % Check if the Data Input box was cnaceled
    return                                                                  % End script if there isn't an input for it to use
end



%% function Data = SonarDataStream2Mat(file)

fileID = fopen(file);
bytes = fread(fileID);
fclose(fileID);

bytes = uint8(bytes);

[~, fileName,ext] = fileparts(file);


%% Create Data Structure
Data = struct('timeStamp', datetime,'position',[], 'utmZone', '', 'heading',[], 'depth',[], 'speed',[], 'portSonar',struct, 'starSonar', struct, 'header','');

Data.header = { strcat('Sonar File:', {' '},fileName, ext);
    'timeStamp: dd/MM/yyyy hh:mm.SS.ssss';
    "position:  Vehicle's position in utms [ Easting, Northing, Altetude ]";
    'utmZone:   Longitude and Latitude bands';
    "heading:   Vehicle's heading [degrees]";
    'depth:     Water depth [meters]';
    "speed:     Vehicle's speed [meters per second]";
    'portSonar: Sonar data structure';
    'starSonar: Sonar data structure';
    "Sonar data structure --> contrast: Sonar's 'Contrast' control setting  in decibels";
    "Sonar data structure --> offset:   Sonar's 'Offset'   control setting  in decibels ";
    "Sonar data structure --> gain:     Sonar's 'Gain'     control setting  in decibels";
    "Sonar data structure --> range:    Sonar's 'Range'    control setting  in Meters";
    "Sonar data structure --> signal:   Signal intensities of acoustic data in decibles ";
    "Sonar data structure --> EchoStrength: echo strength for visulization. Nomalized value, may contain NaN values";
    "Sonar data structure --> SamplePeriod: Sample period for the individual acoustic pings in seconds";
    "Sonar data structure --> Resolution:   Range Resolution for each sample in meters"};


%% Create variable for the incoming data and pre alocate some space for them
preSize = 10000;

heading(preSize)      = 0;        % Prealocate Memeory
position(preSize,1:3) = [0, 0, 0];
utmZone(preSize)      = " ";
speed(preSize)        = 0;
sat(preSize)          = 0;
depth(preSize)        = 0;
dateTime(preSize)     = datetime('now','Format','MMMM d,yyyy HH:mm:ss.SSS');
portSonar(preSize)    = struct('signal',[],'range',[],'speedOfSound',[]);
starSonar(preSize)    = struct('signal',[],'range',[],'speedOfSound',[]);
hardware(preSize)     = struct('Rev',[],'Temp',[],'VIN',[],'IIN',[],'Channel1',[],'Channel2',[],'Channel3',[],'Channel4',[],'Warnings',[],'Alarms',[]);


%% Cycle through the rest of the the data file and extract info
index = 1;
pkgStart = 1;
dateTime_Ref = [];

while pkgStart < length(bytes)
    
    %  Reader Pakage Header -------------------------------------------------------------------
    payloadSize = typecast(bytes(pkgStart  : pkgStart+ 3), 'uint32');       % The size of the entire data packet, not including this header.
    timeStamp   = typecast(bytes(pkgStart+4: pkgStart+ 7), 'uint32');       % The number of milliseconds since boot.
    payloadType = typecast(bytes(pkgStart+8: pkgStart+11), 'uint32');       % This describes the contents of the data packet and how the header will be used.
    miss        = bytes(pkgStart+12 : pkgStart+14);                         % Used for miscellaneous storage. Mainaly the sync package
    chechecksum = bytes(pkgStart+15);                                       % This is a simple 8 bit checksum
    
    % Do checksum -------------------------------------------------------------------------------------------------------------------
    check = xor(bytes(1), bytes(2));
    for ii = 3:16
        check = xor(check, bytes(ii));
    end
    
    if check
        disp('!! ++++++  I no checksum --------- !!!')
    end
    
    
    % Read Payload ---------------------------------------------------------------------------
    payload = bytes( pkgStart+16 : pkgStart +15 + payloadSize);             % Data Package
    
    pkgStart = pkgStart +16 + payloadSize;                                  % Calculate the starting byte of the next data packet
    
    % Interpret payload and assign to data structure -------------------------------------------------------------------------------------------------------
    switch dec2hex(payloadType)
        
        case '534E5200' % Sonar Pakage
            disp('Sonar Pakage')
            
            
        case '53594E43'  % Sync Pakage
            % disp('Sync Pakage')
            if ~isequal(dec2hex(miss), ['AA';'AA';'AA'])                    % Check that the sync packet is correct
                disp('!!!!!!!!  Parsing is Not Syncronized   !!!!!!!!')     % Terminate function if the package does not sync
                return
            end
            dateTime_Ref = datetime(typecast(payload(1:4), 'uint32'),'ConvertFrom','posixtime','Format','MMMM d,yyyy HH:mm:ss.SSS'); % Date and time of the data packate
            % syncInterval = milliseconds(typecast(payload(5:6), 'uint16'));  % Sync packet repetition interval
            lastDatetime = timeStamp;
            
            
        case '4F524E54' % Orientation Pakage
            disp('Orientation Pakage')
            
            
        case '4E415600' % Navigation Pakage
            disp('Navigation Pakage')
            
            
        case '46415400' % Fathometer Pakage
            disp('Fathometer Pakage')
            
            
        case '4D414700' % Magnetometer Pakage
            disp('Magnetometer Pakage')
            
            
        case '4D41524B' % Mark Pakage
            disp('Mark Pakage')
            
            
        case '4E4D4541' % NMEA Pakage
            disp('NMEA Pakage')
            
            
        case '32524E52'
            % disp('Sonar Data V2')
            [portSonar(index), starSonar(index)] = Sonar_Data(payload);
            
            
        case '43554245' % CUBE Data Pakage
            % disp('CUBE Data Pakage')
            hardware(index) = CUBE(payload);
            
            if isempty(dateTime_Ref)
                dateTime(index) = timeStamp;
            else
                dateTime(index) = dateTime_Ref + milliseconds( timeStamp - lastDatetime );
            end
            
            index = index+1;
            
        otherwise
            disp('  ! Unrecognized Payload Type:')
            disp( dec2hex(payloadType) )
            
    end
end


heading(index:preSize)    = [];  % Clean up un ussed space
position(index:preSize,:) = [];
utmZone(index:preSize)    = [];
speed(index:preSize)      = [];
sat(index:preSize)        = [];
depth(index:preSize)      = [];
dateTime(index:preSize)   = [];
portSonar(index:preSize)  = [];
starSonar(index:preSize)  = [];
hardware(index:preSize)   = [];


% Data.heading   = heading;
% Data.position  = position;
% Data.utmZone   = char(utmZone');
% Data.speed     = speed;
% Data.satellite = sat;
% Data.depth     = depth;
Data.timeStamp = datetime( dateTime, 'Format', 'MMMM d,yyyy HH:mm:ss.SSS');
Data.portSonar = portSonar;
Data.starSonar = starSonar;
Data.hardware  = hardware;

disp('       Sonar File Imported')





%% Functions ==================================================================================================================
% Channel specific information -----------------------------------------------------------------------------------------------
function C = CUBE(bytes)

C.Rev  = bytes(1);                              % Stores the revision number for the packet.

C.Temp = typecast( bytes( 2: 5), 'single');                      % This is the system temperature in degrees Celsius.
C.VIN  = typecast( bytes( 6: 9), 'single');                      % This value is the input voltage in Volts DC.
C.IIN  = typecast( bytes(10:13), 'single');                      % This value is the input current in Amps.

C.Channel1 = CUBE_CHANNEL( bytes(14: 35) );                     % CUBE_CHANNEL structures which contain additional channel specific information.
C.Channel2 = CUBE_CHANNEL( bytes(36: 57) );
C.Channel3 = CUBE_CHANNEL( bytes(58: 79) );
C.Channel4 = CUBE_CHANNEL( bytes(80:101) );

C.Warnings = Warn_Alarm( bytes(102) );                          % This field contains the Scout Lite's warning flags.
C.Alarms   = Warn_Alarm( bytes(103) );                          % This field contains the Scout Lite's alarm flags.

end

function S = CUBE_CHANNEL(bytes)

S.Mode        = typecast( bytes( 1: 2), 'uint16');  % This field stores the channel's ping mode (Disabled, Chirp UP, Chirp Down, or CW).
S.CenterFreq  = typecast( bytes( 3: 6), 'single');  % This field stores the center frequency for this channel.
S.Bandwidth   = typecast( bytes( 7:10), 'single');  % This field stores the channel's bandwidth in kHz.
S.PulseLength = typecast( bytes(11:14), 'single');  % This field stores the channel's configured pulse length in microseconds.
S.Duty        = typecast( bytes(15:18), 'single');  % This field stores the channel's transmit power (0.0 to 1.0)
S.Cycles      = typecast( bytes(19:22), 'int32');   % The number of cycles in the sonar pulse at the center frequency if the CW mode is configured for this channel.

end


% Hardware Warnings ----------------------------------------------------------------------------------------------------------
function warn = Warn_Alarm(byte)

switch dec2hex(byte)
    case '0'
        warn = 'NON';
        
    case '1'
        warn = 'Over Temp';
        
    case '2'
        warn = 'Over Volt';
        
    case '3'
        warn = 'Over Current';
        
    otherwise
        disp('Unrecognize Warning')
        disp(dec2hex(byte))
        
end

end


% Sonar Ping -----------------------------------------------------------------------------------------------------------------
function [PortSS, StarSS] = Sonar_Data(bytes)
% Ping data --------------------------------------------------------------------------------------
ChannelCount = bytes(1);                            % Contains a count of the number of sonar channels contained within this data packet.
% pingNumber   = typecast( bytes( 2: 5),'uint32');  % This is the ping number for this sonar. The value increments with every ping (typically every sonar data packet).  This can be used to detect skipped pings etc.
speedOfSound = typecast( bytes( 6: 9),'single');    % Speed of sound in water in meters per second. Typically 1500 m/s.

% reserved      = bytes(10:17);                     % 7 bytes reserved for later use

% Chanel data ------------------------------------------------------------------------------------
sonarType(ChannelCount)    = NaN; % Pre-alocate some memory
sonarID(ChannelCount)      = NaN;
freqHz(ChannelCount)       = NaN;
rangeMs(ChannelCount)      = NaN;
rangeDelayMs(ChannelCount) = NaN;
dataFlags(ChannelCount)    = NaN;
dataType(ChannelCount)     = NaN;
numSamples(ChannelCount)   = NaN;


begin = 17;
for ii = 1:ChannelCount
    sonarType(ii)    = typecast( bytes(begin   : begin+ 1),'uint16');       % Stores the type of sonar used to the collect the data.
    sonarID(ii)      = typecast( bytes(begin+ 2: begin+ 3),'uint16');       % Stores an identifier such as serial number or transducer type for this channel.
    freqHz(ii)       = typecast( bytes(begin+ 4: begin+ 7),'single');       % The center frequency of this channel
    rangeMs(ii)      = typecast( bytes(begin+ 8: begin+11),'single');       % Stores the 1 way slant range in decimal milliseconds.
    rangeDelayMs(ii) = typecast( bytes(begin+12: begin+15),'single');       % Stores the 1 way range delay in decimal seconds. This is the delay from the start of the transmit pulse to the first sample.
    dataFlags(ii)    = typecast( bytes(begin+16: begin+17),'uint16');       % This field stores additional information about this channel.
    dataType(ii)     = typecast( bytes(begin+18: begin+19),'uint16');       % Stores the Sample data type.
    numSamples(ii)   = typecast( bytes(begin+20: begin+21),'uint16');       % The number of bytes for the samples recorded for this channel.
     
    begin = begin + 22;
    
end


% Read Sonar Data ------------------------------------------------------------------------------------------
numSampleBytes = numSamples * 4;        % Get sample size in bytes
sonarType = dec2hex( sonarType );       % Convert sonar type into hex
% dataType  = dec2hex( dataType );        % Convert sonar type into hex
range = rangeMs./1000 * speedOfSound; % Calculate range as a distance


samples(ChannelCount) = NaN; % Pre-alocate some memory

for ii = 1:ChannelCount
    
    last         = begin+numSampleBytes(ii)-1;    
    samplesBytes = reshape( bytes( begin : last ), [numSamples(ii),4]);
    begin        = begin+numSampleBytes(ii);
    
    samples(numSamples(ii)) = NaN;
    
    for jj = 1: numSamples(ii)
        samples(jj) = typecast( samplesBytes(jj,:),'uint32');  % Sonar Samples
    end
    
    switch sonarType(ii,:)
        case '01' % Port SS Daata
            PortSS.signal = samples;
            PortSS.range  = range(ii);
            PortSS.speedOfSound = speedOfSound;
 %           PortSS.EchoStrength = [] ; 
            
        case '02' % Star SS Data
            StarSS.signal = samples;
            StarSS.range  = range(ii);
            StarSS.speedOfSound = speedOfSound;
            
        case '10' % Port SS Data
            PortSS.signal = samples;
            PortSS.range  = range(ii);
            PortSS.speedOfSound = speedOfSound;
            
        case '20' % Star SS Data
            StarSS.signal = samples;
            StarSS.range  = range(ii);
            StarSS.speedOfSound = speedOfSound;
            
        otherwise
            disp('Undefined Sonar Tye')
            disp(sonarType(ii,:))
    end
end

end




