%% Isolate area areound a waypoint

close all
% clear all
clc


%% Load Data
[filepath,name,~] = fileparts( mfilename('fullpath') );                 % Get the file path of the folder
load(fullfile(filepath,'LNH_OutletTower_Watpoint.mat'));                % Load Waypoint
load(fullfile(filepath,'LNH_Region_15.mat'));                           % Load Data
[geoImage, geoData] = geotiffread(fullfile(filepath,'LNH_Damn.tiff'));  % Load geotiff

if ~isa(geoData,'map.rasterref.GeographicCellsReference')               % Check that geoData is in the right format for 'geoshow'
    
    geoData = georasterref('RasterSize',size(geoImage), ...             % Reformat geoData for 'geoshow'
        'RasterInterpretation', geoData.RasterInterpretation, ...
        'LongitudeLimits',geoData.XWorldLimits, ...
        'LatitudeLimits',geoData.YWorldLimits, ...
        'ColumnsStartFrom',geoData.ColumnsStartFrom,...
        'RowsStartFrom',geoData.RowsStartFrom);
end


%% Plot Area data
Figure = figure('Name','Outlet Tower Bathymetry Mapper',...
    'units','normalized',...
    'Position',[0.1  0.1   0.5  0.5],...
    'Color',[1 1 1],...
    'NumberTitle','off',...
    'ToolBar','none' );


hold on
geoshow(geoImage,geoData);

p1 = plot(OutletTower(1),OutletTower(2),'Color','g','Marker','+');

title('Area of Interest')
xlabel('Easting [deg]')
ylabel('Norting [deg]')

legend(p1,'Outlet tower')


%% Get window Size
disp('Enter Window Size in Meters')
windowSize = inputdlg('Enter Window Size [Meters]', 'Input Specs...',[1 20], {'200'});   % Input Diolog box

if isempty(windowSize)
    return
else
    windowSize = str2double(windowSize{1})/2;                                                              % Divide in half, window is created by adding and subtracting this length to the tower's waypoint
end


%% Calculate area around the tower
% Get waypoint in UTMs
[x,y,utmzone] = deg2utm(OutletTower(2),OutletTower(1));

box = [ x + windowSize, y + windowSize;
    x + windowSize, y - windowSize;
    x - windowSize, y - windowSize;
    x - windowSize, y + windowSize;
    x + windowSize, y + windowSize ];

% Convert back to degrees
utmzone = repmat(utmzone,5,1);                                  % Make "utmzone" the right size for "utm2deg"
[window(:,2),window(:,1)] = utm2deg(box(:,1),box(:,2),utmzone); % Arrange Window as [Lon, Lat]


%% Plot Window
p2 = plot(window(:,1), window(:,2));
p2.Color = 'r';
p2.LineWidth = 1;

title('Area of Interest')
xlabel('Easting [deg]')
ylabel('Norting [deg]')
legend([p1,p2],{'Outlet tower','Bathyemetry Area'})

hold off


%% Extract data from the data structure
depths = cat(1,LNH_regionData);
longitude = depths(:,1);
latitude  = depths(:,2);
Data  = depths(:,3);

% Find the data point that are inside the Window
in = inpolygon(longitude, latitude, window(:,1), window(:,2));
longitude = longitude(in);
latitude  = latitude(in);
Data      = Data(in);


%% EKF Filtering
disp('Enter Filter Specs.')
fields = {'Grid cell width','Floor of the map','Window size of EKF','Sigmas'};
answers = inputdlg(fields, 'Input Specs...',[1 20; 1 20; 1 20; 1 20], {'1.0',' -60', '40', '2'} ); % Dialog box for user to enter EKF filter specs

if isempty(answers)     % Use Default answers if diolog box is canceled
    res        =   1.0; % Grid cell width
    baseFloor  = -60;   % Floor of the map --> This will defalt to the minimum data value if the argument is not passed into the Bathymetry function.
    filterSize =  40;   % Window size of EKF
    sig        =   2;   % Sigma threshold for low pass filter
else
    res        = str2double(answers{1}); % Grid cell width
    baseFloor  = str2double(answers{2}); % Floor of the map --> This will defalt to the minimum data value if the argument is not passed into the Bathymetry function.
    filterSize = str2double(answers{3}); % Window size of EKF
    sig        = str2double(answers{4}); % Sigma threshold for low pass filter
end

% Run EKF to create the bathymetry map
[Bathymetry, LON,  LAT, Var] = EKF_2D( longitude, latitude, Data,'GridCell', res,'Window', filterSize, 'sigma', sig); % ,'Floor', baseFloor );


%% Convert coordinates to Meters to scale axes equaly
%   * EKF_2D is set up to work in degrees

xOrg = min(LON(:));     % find Origin of map
yOrg = min(LAT(:));

xdim = LON(1,:)';       % Get map increments
ydim = LAT(:,1);

xOrg = repmat(xOrg,length(LON),1);  % Make "xOrg" the right size for "vdist"
yOrg = repmat(yOrg,length(LON),1);

[rows, colums] = size(LON);         % Get the size of the Longitude an latitude mesh to mach the new x-t intervals

xInt = vdist(yOrg(1:colums), xOrg(1:colums), yOrg(1:colums), xdim);         % Calculate the distances of the LAT LON grids in meters
yInt = vdist(yOrg(1:rows),   xOrg(1:rows),   ydim,           xOrg(1:rows));

[xGrid, yGrid] = meshgrid(xInt, yInt);      % Make new goe-referencing meshies in meters

xGrid = xGrid - windowSize;                 % Shift the new grids to be centered on the outlet tower
yGrid = yGrid - windowSize;


%% Plot the Bathymetry Map
figure('Name','Bathymetry Model', ...
    'Color',[1 1 1],...
    'NumberTitle','off' );

hold on
axis equal

%-- Create Surface from the Bathymetry Data -------------------------------
surface(xGrid, yGrid, Bathymetry,'edgecolor', 'none');
colormap(copper);
caxis([(min(Bathymetry(:)-10)), 0])
l = light('Position',[0 0 3000],'Style','infinite');
lighting gouraud
material dull

%-- Overlay Contours ------------------------------------------------------
v = 0:-2.5:floor(min(Bathymetry(:)));
[c,h] = contour3(xGrid, yGrid,Bathymetry,v, 'k');
v = 0:-10:floor(min(Bathymetry(:)));
clabel(c,h,v,'FontSize',8)
colorbar

% Plot location of Outlet Tower -------------------------------------------
z = Bathymetry(round(rows/2), round(colums/2));
s = scatter3(0,0,z,36,'g');
hold off

title('Interpolated Batymetry')
xlabel('Easting [meters]')
ylabel('Norting [meters]')
zlabel('Depth [meters]')
legend(s,'Location of Outlet Tower')

view(-80,43)

disp('Done')



%% Functions ==========================================================================================
function AreaWindow(varargin)

Figure = figure('Name','Select Area', ...                               % Figure windo
    'units','normalized',...
    'Position',[0.5  0.5   0.3  0.3],...
    'NumberTitle','off',...
    'ToolBar','none', ...
    'MenuBar','none' );

uicontrol(Figure,'Style','text','units','normalized', ...           % Labe --> Enter waypoint
    'Position', [0.01  0.9  0.9 0.05], ...
    'HorizontalAlignment','left',...
    'FontSize', 12, ...
    'BackgroundColor',[ 0.9294    0.6902    0.1294 ], ...
    'String','Enter a Waypoint in the Center of the Isolated Area or Load Waypoint from .mat variable');

uicontrol(Figure,'Style','edit','units','normalized', ...           % Enter Longitude
    'Position', [0.1  0.8  0.5 0.07], ...
    'Tag','Lon', ...
    'String', 'Longitude', ...
    'Callback', @GetArea );

uicontrol(Figure,'Style','edit','units','normalized', ...           % Enter Latitude
    'Position', [0.1  0.7  0.5 0.07], ...
    'Tag','Lat', ...
    'String', 'Latatude', ...
    'Callback', @GetArea );


uicontrol(Figure,'Style','pushbutton', 'units','normalized',...     % Load Waypoint from .mat
    'Position', [0.65 0.8 0.3 0.07], ...
    'Tag','Mat', ...
    'String', 'Load Waypoint from .mat', ...
    'Callback', @GetArea );

uicontrol(Figure,'Style','text','units','normalized', ...           % Label -> Window Size
    'Position', [0.01  0.5  0.6 0.05], ...
    'BackgroundColor',[ 0.3020    0.7490    0.9294 ], ...
    'HorizontalAlignment','left',...
    'FontSize', 12, ...
    'String','How Big should the area be?');

uicontrol(Figure,'Style','edit','units','normalized', ...           % Enter window Size
    'Position', [0.1  0.4  0.5 0.07], ...
    'Tag', 'Window', ...
    'String', 'Window Size', ...
    'Callback', @GetArea );

deleteFigure = @(a,b) delete(Figure);                               % calbak function to delete figure window

uicontrol(Figure,'Style','pushbutton', 'units','normalized',...     % Done button
    'Position', [0.7 0.1 0.2 0.07], ...
    'String', 'Done', ...
    'Callback', deleteFigure );

waitfor(Figure)
disp('Window Closed')

end


function inPoints = IsolateArea(varargin)
global Area;


inPoints = nan;


[x,y,utmzone] = deg2utm(Area{1},Area{2});

x
y

windowSize = Area{3}

box = [ x + windowSize, y + windowSize;
    x + windowSize, y - windowSize;
    x - windowSize, y - windowSize;
    x - windowSize, y + windowSize;
    x + windowSize, y + windowSize ];

% Convert back to degrees
utmzone = repmat(utmzone,5,1);                                  % Make "utmzone" the right size for "utm2deg"
[window(:,2),window(:,1)] = utm2deg(box(:,1),box(:,2),utmzone); % Arrange Window as [Lon, Lat]

window

hold on
p1 = plot(Area{1},Area{2},'Color','g','Marker','+');
p2 = plot(window(:,1), window(:,2));
p2.Color = 'r';
p2.LineWidth = 1;

hold off

title('Area of Interest')
xlabel('Easting [deg]')
ylabel('Norting [deg]')
legend([p1,p2],{'Outlet tower','Bathyemetry Area'})

end


%% Call Back Functions
function GetArea(varargin)
global UserInputs
global Area



switch varargin{1}.Tag
    
    case 'Lon'
        Area{1} = str2double( varargin{2}.Source.String )
        
        
    case 'Lat'
        Area{2} = str2double( varargin{2}.Source.String )
        
        
    case 'Mat'
        try
            [DataFileName,DataFilePath] = uigetfile('*.mat','Select MATLAB Data File',UserInputs.Preferences('WayPoint'));
        catch
            [DataFileName,DataFilePath] = uigetfile('*.mat','Select MATLAB Data File');
        end
        
        if DataFileName == 0                                % Check if the Data Input box was cnaceled
            disp('Input box Cancelled')
            return                                          % End script if there isn't an input for it to use
            
        else
            dataFile = fullfile(DataFilePath,DataFileName); % Get full file path
            WayPoint = load(dataFile);                % Load Waypoint
            fnames = fieldnames(WayPoint);
            
            Area{1} = WayPoint.(fnames{1})(1);
            Area{2} = WayPoint.(fnames{1})(2);
            
            UserInputs.Preferences('WayPoint') = dataFile;  % Save file path
            UserInputs.SavePreferences();
        end
        
        
    case 'Window'
        Area{3} = str2double( varargin{2}.Source.String )
        
        
    otherwise
        varargin{1}.Tag
        
end

end


