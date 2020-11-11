% Jacob Anderson
% 12-05-2017
% jaanderson@fortlewis.edu

% This script is called by AkaboticsGUI.m
% This scrip will create and modify an Occupangy Grid that covers a geografical area of interest.
% The Occupancy is to be used for environmentla niche modeling and robotic navigation.

function MakeOcupancyGrid(fromPanel)
global MainWindow;
global Model;                               % Area model being developed by this app
global geoImage

global shape;               % List of waypoints provided by the user. Waypoints outline the area of interest
global LON;                 % Longitude grid for geo-referancing
global LAT;                 % Latitude grid for geo-referancing


try
    BW = evalin('base', 'BW');  % See if a blak and white maks is available from color threasholding
catch
    BW = [];                    % Create an empty array if it isn't there
end

% Create the ocupancy Grid-------------------------------------------------------------------------------
if ~isempty(shape) && isempty(LON) %----------------------------> A shape file exists but there is no geotiff
    x = linspace( min(shape(:,1)), max(shape(:,1)), 500 );      % Get geographic limits from shape file
    y = linspace( min(shape(:,2)), max(shape(:,2)), 500 );
    [LON, LAT] = meshgrid(x,y);                                 % Create Latitude and Longitude grids for goe-refrencing
    
    BW = inpolygon(LON,LAT,shape(:,1),shape(:,2));              % Create a mask of the area enclosed in the shape file
    OcupancyGrid = BW;                                          % Coppy the mask to the occupancy Grid.
    
elseif ~isempty(shape) && ~isempty(LON) %-----------------------> Shape file and Geotif exists
    BW = inpolygon(LON,LAT,shape(:,1),shape(:,2));              % Create a mask of the area enclosed in the shape file
    OcupancyGrid = BW;                                          % Coppy the mask to the occupancy Grid.
    
elseif isempty(shape) && ~isempty(geoImage) %-------------------------> There is a color threshold from the Geotiff and no Shanpe file
    
    lab_he = rgb2lab(geoImage);
    ab = lab_he(:,:,2:3);
    ab = im2single(ab);
    nColors = 3;
    pixel_labels = imsegkmeans(ab,nColors,'NumAttempts',3);
    BW = pixel_labels==2;
    BW = ~BW;
    
    h = fspecial('average',[10 10]);                            % Create an averaging filter to smoth edges
    image = imfilter(BW,h);                                     % Applie the filter
    
    OcupancyGrid = bwareaopen(image,200);                       % Fill in the out of bounds areas
    OcupancyGrid = bwareaopen(~OcupancyGrid,150);               % Fill in the in bounds areas
    OcupancyGrid = ~OcupancyGrid;                               % Re-invert the ocupancy grid so that the in bounds areas are represented by 1s
    
else %----------------------------------------------------------> Some how we got here without the requiered data. 
    disp('Somthings Missing')                                   % Disply error message
    ST = dbstack(1);
    disp(ST(1).name);
    return                                                      % Terminate the script
end

[Boundaries, Regions, ~, parentChild] = bwboundaries(OcupancyGrid,'holes');   % Identify Regions and Boundaries of the Black areas

Model(1).OcupancyGrid = OcupancyGrid;
Model(1).Boundaries   = Boundaries;
Model(1).Regions      = Regions;
Model(1).parentChild  = parentChild;
Model(1).LAT          = LAT;
Model(1).LON          = LON;


MainWindow.BackGround('LON',LON,'LAT',LAT);                     % Send the Longitude and Latatude mesh to the GUI
MainWindow.Addlayer('OcupancyGrid', OcupancyGrid, fromPanel);   % Dispaly OcupancyGrid on the GUI window

MainWindow.Header.String = 'Refine Occupancy Grid by Adding and Removing Area or Adjusting the Boundary';   % Prompt to select data

% Create a button panel on the GUI to organize ocupancy grid tools ===============================================================================================
panalName = 'Occ. grid';
MainWindow.GUIpanels(panalName,'make','position',[0.7 0 0.1 1],'backgroundcolor',[0.9294  0.6902  0.1294]);  % Create putton panel for primary GUI functions

% Make buttons in the button panel {Buttin type,  Button Name, Callback function and inputs, Callbacl function for slider listeners
MainWindow.GUIbuttons( panalName,'add', ...
    {'pushbutton','Adjust Grid',  {@AddjustBoundaries, 'addjust'}, [], 'Addjusts the boundary of the Grid\nMove the nodes around\nDelete unwanted nodes\nDoulbe click inside boundary to accept';
     'pushbutton', 'Add Area',     {@AddRegion,         'add'},     [], 'Add area to the Grid\nLeft click to close the polygon\nDouble click inside polygon to accept';
     'pushbutton', 'Remove Area',  {@AddRegion,         'remove'},  [], 'Remove area from the Grid\nLeft click to close the polygon\nDouble click inside polygon to accept';
     'pushbutton', 'Show Regions', {@ShowRegions,       panalName}, [], '';
     'pushbutton', 'Done',         {@LoadOccGric,       panalName}, [], '';
    });

end


%% Callback Functions
% Addjust boundaries -----------------------------------------------------------------------------------------------------------
% Create an interactive polygon around the Occupancy Grid for the user to maniplulate its size and shape
function AddjustBoundaries(~,~,command)
% Varargin: {'UIControl', 'matlab.ui.eventdata.ActionData',command}
global MainWindow;
global Model
global LAT
global LON

stats = regionprops(Model(1).Regions,'all');                                   % Open all regionprops options for the regions identified
area = cat(1,stats.Area);                                                   % Find the boundary that corresponds to the Occupancy grid
[~,ii] = max(area);

pPoints = Model(1).Boundaries{ii};

[latout,lonout] = reducem( LAT(pPoints(:,1),1), LON(1,pPoints(:,2))',0.00006);  % Boundaries will contain way tooo many referance points for the interactive polygon. This reduces the number points

pos = MainWindow.InteractiveLayer( lonout, latout);                         % Create interactive polygon in the GUI and wait for it to be closes
                                                                            % pose is a list of waypoints from the vertecies of the polygon
NewBoundary(pos,command);                                                   % Process the new region

end


% AddRegion ----------------------------------------------------------------------------
% Get a region to add or subtract from the ocupancy grid via an interactive poly gon
function AddRegion(~,~,command)
% Varargin: {'UIControl', 'matlab.ui.eventdata.ActionData', 'panalName',command}
global MainWindow;
pos = MainWindow.InteractiveLayer;            % Create interactive polygon in the GUI and wait for it to be closes
NewBoundary(pos,command);                     % Process the new region
end


% Display the regions from the blak and with masking --> This is more for the desiger's use and may be removed in the future
function ShowRegions(~,~,panelName)
global MainWindow;
global Model;
MainWindow.Addlayer('Regions', Model(1).Regions, panelName);
end


% Load Ocupancy Grid -----------------------------------------------------------------
% The user is finished creating the occupancy grid and ready to move on
function LoadOccGric(~,~,panelName)
% Varargin: {'UIControl', 'matlab.ui.eventdata.ActionData', 'panalName' }
global MainWindow;
global Model;
global DATA

stats = regionprops(Model(1).Regions,'all');                % Open all regionprops options for the regions identified
area = cat(1,stats.Area);                                   % Find the boundary that corresponds to the Occupancy grid
[~,ii] = max(area);

Boundaries = Model(1).Boundaries{ii};

[Lon, Lat] = XYtoLatLon(Boundaries(:,2),Boundaries(:,1));

Model(1).Boundaries = [Lon, Lat];
MainWindow.Addlayer('Boundaries',[Lon, Lat], panelName);    % Dispaly OcupancyGrid on the GUI window

MainWindow.GUIpanels(panelName,'delete');                   % Close the button panel for the occupancy grid tools

MainWindow.Header.String = 'Select Data Input';             % Prompt to select data


if ~isempty(DATA)
    answer = questdlg('Apply Boundary to Data Points?', ...
        'Apply Boundary to Data Points?', ...
        'Yes','No','No');
    
    switch answer
        case 'Yes'
%             DATA = GoeFilterDATA(DATA, Boundaries);
%             
%             dataPoints = [];
%             
%             for ii = 1: length(DATA)
%                 
%                 dataPoints = [dataPoints; cat(1,DATA{ii}.vehicle)];                                               % Extrace data waypointds from the data structure for visulizations
%                 
%             end
%             MainWindow.Layers('DataPoints') = dataPoints;
%             MainWindow.Plot;
            
        otherwise
    end
end
% Save Boundary
% answer = questdlg('Export Boundary?', ...
% 	'Expor boundary', ...
% 	'Yes','No','No');
% 
% switch answer
%     case 'Yes'
%          uisave('Boundaries',[num2str(yyyymmdd(datetime('today'))),'_Boundaries']);
%         
%     otherwise
% end

end



%% Boundary Work ------------------------------------------------------------------
% Impliment the changes from the interactive polygons
function NewBoundary(pos,command)
% Varargin: {'UIControl', 'matlab.ui.eventdata.ActionData', 'panalName',command}
global MainWindow;
global Model;

% Creat Mask of the area inside the polygon
[m,n] =size(Model(1).OcupancyGrid);
[row,column] = getIndices(pos, m, n);  % Get row and column indecies of the polygon vertecies
BW = poly2mask(column, row, m, n);     % Creat Mask of the area inside the polygon


% Combine the mask with the ocupancy grid
switch command
    
    case 'addjust'
        
        stats = regionprops(Model(1).Regions,'all');       % Open all regionprops options for the regions identified
        
        area = cat(1,stats.Area);
        [~,ii] = max(area);
        
        enclosed_boundaries = find(Model(1).parentChild(:,ii));
        if ~isempty(enclosed_boundaries)
            for i = 1:length(enclosed_boundaries)
                Points = Model(1).Boundaries{enclosed_boundaries(i)};
                BW2 = poly2mask(Points(1), Points(2), m, n);     % Creat Mask of the area inside the polygon
                BW(BW2 == true) = false;
            end
        end
        
        Model(1).OcupancyGrid = BW;
        
    case 'add'
        Model(1).OcupancyGrid = Model(1).OcupancyGrid | BW;
            
    case 'remove'
        Model(1).OcupancyGrid(BW == true) = false;
        
    otherwise
        error('invalid input')
        
end

MainWindow.Layers('OcupancyGrid') = Model(1).OcupancyGrid;
MainWindow.Plot;

[Boundaries, Regions, ~, parentChild] = bwboundaries(Model(1).OcupancyGrid,'holes');   % Identify Boundaries of the Black areas

Model(1).Boundaries  = Boundaries;
Model(1).Regions     = Regions;
Model(1).parentChild = parentChild;

MainWindow.Header.String = 'Refine Occupancy Grid by Adding and Removing Area or Adjusting the Boundary';                             % Prompt to select data

end




%% Supporitng Functions ==================================================================

% Determine the row,column indicies of a cell given an x,y location
function [row,column] = getIndices(L, numRows, numCols)
global LAT;
global LON;

[m,~] = size(L);                     % Deterimin the number of enteries in L for vectorization

lat0(1:m,1)  = min(LAT(:));          % Get the corners of the map
long0(1:m,1) = min(LON(:));
latD(1:m,1)  = max(LAT(:));
longD(1:m,1) = max(LON(:));

lat  = L(:,2);                       % Extract lat and lon from input
long = L(:,1);

x = vdist(lat0,long,lat0,long0);     % from longitude to x in meters
y = vdist(lat,long0,lat0,long0);     % from latitude  to y in meters

x(isnan(x)) = 0;                     % Filter out NaNs
y(isnan(y)) = 0;

xMax = vdist(lat0(1),longD(1),lat0(1),long0(1));    % Get length of the map, x dimention must be inverted becaus of the way imaging software plots
yMax = vdist(latD(1),long0(1),lat0(1),long0(1));

row    = numRows - floor(y ./( yMax / numRows));    % Y indecies
column = floor(x ./( xMax / numCols));              % X indecies

row(row<=0) = 1;                                    % Make sure there are no 0 indecies
column(column <=0) = 1;
end


% Converts latitude to an x value, and longitude to a y value
%       L is the 2x1 vector with latitude, longitude
%       L0 is the 2x1 vector with the latitude, longitude of the coordinate frame origin
% function [x,y] = getXYCoordinates(L)
% global LAT;
% global LON;
% 
% 
% lat  = L(:,2);
% long = L(:,1);
% 
% [m,~] = size(L);
% 
% lat0(1:m,1) = min(LAT(:));
% long0(1:m,1) = min(LON(:));
% 
% x = vdist(lat0,long,lat0,long0); % from longitude to x in meters
% y = vdist(lat,long0,lat0,long0); % from latitude  to y in meters
% 
% x(isnan(x)) = 0;
% y(isnan(y)) = 0;
% end


% Get Lon and Lat from XY coordinates
function [Lon, Lat] = XYtoLatLon(X, Y)
global LAT;
global LON;

Lon = LON(1,X)';
Lat = LAT(Y,1);

end

