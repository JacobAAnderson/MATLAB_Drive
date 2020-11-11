%% Jacob Anderson
% Robotic GNOME Laboratory
% Fort Lewis College
% Durango, CO, 81301
%
% 2017
%
% All Rights Reserved.
%
% Contact: Ryan N. Smith - rnsmith@fortlewis.edu
%
% This is open software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version. It would be appreciated to provide
% credit to the originating source when appropriate.
%
%  This is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You can find the GNU General Public License at <http://www.gnu.org/licenses/>.
%
% -----------------------------------------------------------------------------------------------
% EcoMapperLog2Mat reads in a .log file created by an Ecomapper AUV
%
% Inputs:
%   file --> The full file path to the  .log file
%   varargin --> (optional) Filter data points based on DVL fix
%
% Outputs:
%   Data --> A structure that contain the imported data and metadata
% -----------------------------------------------------------------------------------------------



%% EcoMapperLog2Mat  ==========================================================================================================================================================
function Data = EcoMapperLog2Mat(file, filter)

[~, fileName, ext] = fileparts(file);

fprintf(' * Importing Ecomapper .log file %s\n', fileName);

%% Create Data Structure
Data = struct('timeStamp',[],'vehicle',[],'bathymetry',[],'wqData',[], 'attitude', [], 'tide', [], 'filteredOut', 0,'header','','log','' );

Data.header = {'  timeStamp: Date and time [MM.dd.yyyy HH:mm:ss.SS]';
               '    vehicle: Latitude [deg.dd], Longitude [deg.dd], depth from surface [m]';
               '      bathy: Latitude [deg.dd], Longitude [deg.dd], depth [m]';
               '     wqData: Temp [c], SpCond [mS/cm], Sal [ppt], pH, Turbid + NTU, Chl [ug/L], BGA-PC [cells/ml], ODO [%], ODO [mg/L]';
               '   attitude: Roll [deg], Pitch [deg], Yaw [deg], Altitude [m], Speed [m/s]'
               'filteredOut: Number of data points filtered out due to no DVL fix'; 
               '     header: Explains the entries in this data structure';
               '        log: Name of the .log file that the data came from';
               'Known Stats: Altimeter: 0 mean, 0.059687[m] std' };
            

Data.log = [fileName,ext];    % Add name of the log file to the data structure

Data.timeStamp = NaT;         % Default timestamp value incase the file is unreadable

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

% Open the file, scan in the text and then close the file
formatSpec = '%f%f%{HH:mm:ss.SS}D%{MM/dd/yyyy}D%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s%f%f%f%f%f%f%f%f%f%f%f%f%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%{MM/dd/yyyy}D%{HH:mm:ss}D%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
delimiter = ';';
startRow = 2;
try
    fileID = fopen(file,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
    fclose(fileID);
    
catch exception                 % File Cannot be read --> Skip it and return empty struct
    warning off backtrace
    warning('Error Reading log file \n\t\t%s \n\t\t%s \n\t\tError in: %s, Line: %d \n\n\t\tSuspect bad log file, skipping it \n\n',exception.identifier,exception.message,exception.stack(1).name,exception.stack(1).line)
    warning on backtrace
    return
    
end


% Extrace the data from their cell arrays and eliminate the last entry incase the last line of the log file is incomplet
L = length(dataArray{1});

Latitude  = dataArray{1};
Longitude = dataArray{2};


if nargin > 1 && strcmpi(filter,'dvl filter')
    try
        inpoints = logical(998 > dataArray{47});
        
    catch exception
        fprintf('\t--> ERROR: %s \n',exception.message)
        DVLfix = dataArray{47};
        I = isnan(DVLfix);
        DVLfix(I) = 0;
        DVLfix = logical(998 > DVLfix);
        DVLfix(L) = false;
        
        inpoints = DVLfix;
    end
    
    if sum(inpoints) == 0
        fprintf('\t--> This file failed to pass the filters, ommiting all data points \n\n')
        return
        
    elseif sum(inpoints) < L
        fprintf('\t--> %5d out of %5d points filtered out for no DVL fix \n\n',L-sum(inpoints),L)
        Data(1).filteredOut = L-sum(inpoints);                                       % Make a note of how many data points were filtered out
        
    else
        % Do Nothing... Everything is Great :)
    end
    
else
    inpoints = true(L,1);
end

inpoints(L) = 0;                            % Remove last line from the log file incase it is incompleat


Latitude     = Latitude(inpoints);
Longitude    = Longitude(inpoints);

Time         = dataArray{3};        Time         = Time(inpoints);
Date         = dataArray{4};        Date         = Date(inpoints);
timeStamp    = datetime([Date.Year  Date.Month  Date.Day  Time.Hour  Time.Minute Time.Second], 'Format', 'MM.dd.yyyy HH:mm:ss.SS');

% Vehicle Data -----------------------------------------------------------------
Yaw          = dataArray{11};       Yaw          = Yaw(inpoints);
Pitch        = dataArray{12};       Pitch        = Pitch(inpoints);
Roll         = dataArray{13};       Roll         = Roll(inpoints);
DFSDepthm    = dataArray{15};       DFSDepthm    = DFSDepthm(inpoints);
DTBHeightm   = dataArray{16};       DTBHeightm   = DTBHeightm(inpoints);
Speed        = dataArray{28};       Speed        = Speed(inpoints);

% Water Quality Data -----------------------------------------------------------
TempC        = dataArray{57};       TempC        = TempC(inpoints);
SpCondmScm   = dataArray{58};       SpCondmScm   = SpCondmScm(inpoints);
Salppt       = dataArray{59};       Salppt       = Salppt(inpoints);
%  Depthfeet    = dataArray{60};       Depthfeet    = Depthfeet(inpoints);
pH           = dataArray{61};       pH           = pH(inpoints);
%  pHmV         = dataArray{62};       pHmV         = pHmV(inpoints);
TurbidNTU    = dataArray{63};       TurbidNTU    = TurbidNTU(inpoints);
ChlugL       = dataArray{64};       ChlugL       = ChlugL(inpoints);
BGAPCcellsmL = dataArray{65};       BGAPCcellsmL = BGAPCcellsmL(inpoints);
ODOsat       = dataArray{66};       ODOsat       = ODOsat(inpoints);
ODOmgL       = dataArray{67};       ODOmgL       = ODOmgL(inpoints);




% Correct Bathymetry Data for the orientation of the vehicle ---------------------------------------------------------------------------------------------------------------------------
[x,y,utmzone] = deg2utm(Latitude,Longitude);                                % Vehilces position in utm --> [x,y,utmzone] = deg2utm(Lat,Lon)

ecomapper = num2cell([x,y,DFSDepthm]',1);                                   % Vehicles position in 3 space

p_eco = arrayfun( @(z) [0;0;z],                                                            DTBHeightm,           'UniformOutput',false); % Depth reading as a vector under the vehicle
Rx    = arrayfun( @(roll) [ 1 0 0; 0 cosd(roll) -sind(roll); 0 sind(roll) cosd(roll)],     -Roll,                'UniformOutput',false); % X rotation matrix with role angel
Ry    = arrayfun( @(pitch)[ cosd(pitch) 0 sind(pitch); 0 1 0; -sind(pitch) 0 cosd(pitch)], -Pitch,               'UniformOutput',false); % Y rotation matrix with pitch angle
Rz    = arrayfun( @(yaw)  [cosd(yaw) -sind(yaw) 0; sind(yaw) cosd(yaw) 0; 0 0 1],          -Yaw,                 'UniformOutput',false); % Z rotation matrix with compass heading
R     =  cellfun( @(a, b, c) a*b*c,                                                        Rz, Ry, Rx,           'UniformOutput',false); % Compounded rotation matrix
P_w   =  cellfun( @(r,eco,p) [ r, eco; 0 0 0 1] * [p;1],                                   R, ecomapper', p_eco, 'UniformOutput',false); % Coordinate transormation

P_w = cell2mat(P_w);
P_w = reshape(P_w,[4,sum(inpoints)])';

[Lat_ ,Lon_] = utm2deg(P_w(:,1),P_w(:,2),utmzone);                          % Convert the bottom point's coordinates into lat lon


% Convert vehicle speed from knots to meters per second -------------------------------------------------------------------------------------------------------------------------------
Speed = Speed * 0.514444;

% Assign Values to Data Structure --------------------------------------------------------------------------
Data.timeStamp  = timeStamp;                                                % Date and time of the measurment
Data.vehicle    = [Latitude, Longitude, DFSDepthm];
Data.bathymetry = [Lat_, Lon_, P_w(:,3)];
Data.wqData     = [TempC, SpCondmScm, Salppt, pH, TurbidNTU, ChlugL, BGAPCcellsmL, ODOsat, ODOmgL];
Data.attitude   = [Roll, Pitch, Yaw, DTBHeightm, Speed];
Data.tide       = zeros(size(timeStamp));
end




