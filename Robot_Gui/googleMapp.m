clear all
close all
clc


lon = -107.8982160;
lat =   37.2213890;
zoomlevel = 13;
scale = 2;
width = 640;
height = 640;
maptype = 'satellite';

filename = 'tmp.jpg';

[filepath,name,~] = fileparts(mfilename('fullpath'));
filepath = fullfile(filepath, 'Goolge');

if ~exist(filepath,'dir')                              % Create the directory for the saved inputs if it does not exist
    disp('Making directory for saved inputs')
    mkdir(filepath)
end

filepath = fullfile(filepath, filename);

preamble = 'http://maps.googleapis.com/maps/api/staticmap';
location = ['?center=' num2str(lat,10) ',' num2str(lon,10)];
zoomStr = ['&zoom=' num2str(zoomlevel)];
sizeStr = ['&scale=' num2str(scale) '&size=' num2str(width) 'x' num2str(height)];
maptypeStr = ['&maptype=' maptype ];
format = '&format=jpg';
sensor = '&sensor=false';

url = [preamble location zoomStr sizeStr maptypeStr format sensor];

try
    urlwrite(url,filepath);
catch Error% error downloading map
    
    disp(Error)
    
    warning(['Unable to download map form Google Servers.\n' ...
        'Matlab error was: %s\n\n' ...
        'Possible reasons: missing write permissions, no network connection, quota exceeded, or some other error.\n' ...
        'Consider using an API key if quota problems persist.\n\n' ...
        'To debug, try pasting the following URL in your browser, which may result in a more informative error:\n%s'], lasterr, url);
    return
end

%%
imag = imread(filepath);

% Calculate a meshgrid of pixel coordinates in EPSG:900913
width = size(imag,2);
height = size(imag,1);

centerPixelY = round(height/2);
centerPixelX = round(width/2);

[centerX,centerY] = latLonToMeters(lat, lon ); % center coordinates in EPSG:900913

curResolution = initialResolution / 2^zoomlevel / scale / resize; % meters/pixel (EPSG:900913)

xVec = centerX + ((1:width)-centerPixelX) * curResolution; % x vector
yVec = centerY + ((height:-1:1)-centerPixelY) * curResolution; % y vector

[xMesh,yMesh] = meshgrid(xVec,yVec); % construct meshgrid

% convert meshgrid to WGS1984
[lonMesh,latMesh] = metersToLatLon(xMesh,yMesh);

% Next, project the data into a uniform WGS1984 grid
uniHeight = round(height*resize);
uniWidth = round(width*resize);

latVect = linspace(latMesh(1,1),latMesh(end,1),uniHeight);
lonVect = linspace(lonMesh(1,1),lonMesh(1,end),uniWidth);

[uniLonMesh,uniLatMesh] = meshgrid(lonVect,latVect);

uniImag = zeros(uniHeight,uniWidth,3);

% Fast Interpolation to uniform grid
uniImag =  myTurboInterp2(lonMesh,latMesh,imag,uniLonMesh,uniLatMesh);




%% Supporting functions =================================================================================
function [lon,lat] = metersToLatLon(x,y)
% Converts XY point from Spherical Mercator EPSG:900913 to lat/lon in WGS84 Datum
originShift = 2 * pi * 6378137 / 2.0; % 20037508.342789244
lon = (x ./ originShift) * 180;
lat = (y ./ originShift) * 180;
lat = 180 / pi * (2 * atan( exp( lat * pi / 180)) - pi / 2);
end

function [x,y] = latLonToMeters(lat, lon )
% Converts given lat/lon in WGS84 Datum to XY in Spherical Mercator EPSG:900913"
originShift = 2 * pi * 6378137 / 2.0; % 20037508.342789244
x = lon * originShift / 180;
y = log(tan((90 + lat) * pi / 360 )) / (pi / 180);
y = y * originShift / 180;

end

function ZI = myTurboInterp2(X,Y,Z,XI,YI)
% An extremely fast nearest neighbour 2D interpolation, assuming both input
% and output grids consist only of squares, meaning:
% - uniform X for each column
% - uniform Y for each row
XI = XI(1,:);
X = X(1,:);
YI = YI(:,1);
Y = Y(:,1);

xiPos = nan*ones(size(XI));
xLen = length(X);
yiPos = nan*ones(size(YI));
yLen = length(Y);
% find x conversion
xPos = 1;
for idx = 1:length(xiPos)
    if XI(idx) >= X(1) && XI(idx) <= X(end)
        while xPos < xLen && X(xPos+1)<XI(idx)
            xPos = xPos + 1;
        end
        diffs = abs(X(xPos:xPos+1)-XI(idx));
        if diffs(1) < diffs(2)
            xiPos(idx) = xPos;
        else
            xiPos(idx) = xPos + 1;
        end
    end
end
% find y conversion
yPos = 1;
for idx = 1:length(yiPos)
    if YI(idx) <= Y(1) && YI(idx) >= Y(end)
        while yPos < yLen && Y(yPos+1)>YI(idx)
            yPos = yPos + 1;
        end
        diffs = abs(Y(yPos:yPos+1)-YI(idx));
        if diffs(1) < diffs(2)
            yiPos(idx) = yPos;
        else
            yiPos(idx) = yPos + 1;
        end
    end
end
ZI = Z(yiPos,xiPos,:);
end
