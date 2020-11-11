%% Get a map from a Web Map Server
%____________________________________________________________________________________________________________________________________________________
% Default WMS 
% EOX::Maps   --> https://maps.eox.at
% EOX maps are cloudless Sentinal 2 satalite images that have been stiched together
% EOX offers a free Web Map server at: 'https://tiles.maps.eox.at/wms'
%   * Examplez: Basic qurry URL
%     'https://tiles.maps.eox.at/wms?service=wms&request=getmap&version=1.1.1&layers=s2cloudless-2018&bbox=-132.5604992025939,38.27584446206717,-119.42085076509389,51.41549289956717&width=4096&height=4096&srs=epsg:4326'
%____________________________________________________________________________________________________________________________________________________



function [geoImage, geoData] = GetMapfromWMS( lon, lat )


wmsurl = "https://tiles.maps.eox.at/wms?service=wms&request=getmap&version=1.1.1&layers=s2cloudless-2018";  % Base url

bbox   = strcat("&bbox=", ...
                string(min(lon(:))), ',', string(min(lat(:))), ',', ... 
                string(max(lon(:))), ',', string(max(lat(:))) );                                  % Geograpic Image Area
            
width  = "&width=1024";                                                                                     % Image width
height = "&height=768";                                                                                    % Image height
srs    = "&srs=epsg:4326";                                                                                  % ????

url = strcat(wmsurl,bbox,width,height,srs);


disp(' ')
disp('________________________________________________________________________')
disp('Geting Map From:')
disp(url)
disp('________________________________________________________________________')
disp(' ')


try
    [geoImage, geoData] = wmsread( url );                                      % Request map from the server
catch Error
    disp(Error)
    [geoImage, geoData] = wmsread( url );                                      % Request map from the server
end

% Make Sure the image is in the right format
geoImage = im2uint8(geoImage);

return

%% Find a Web Map Server that provides National Agriculture Imagery Program (NAIP) imagery
%   * Other option, not in use
%   * This should return two servers. Choose the first one. Its more recent

% layer = 'naip';                                                                                           % National Agriculture service layer
% wmsSources = wmsfind(layer,'Version','online','Latlim',latLimit,'Lonlim',lonLimit);                       % Find Web Map Servers that provide that layer
% [geoImage, geoData] = wmsread(wmsSources(1),'Latlim',latLimit,'Lonlim',lonLimit,'TimeoutInSeconds',180);  % Request Map from server
