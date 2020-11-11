
% ROMS Resources =====================================================================================================
% Places to get data from:
%   * http://ingria.coas.oregonstate.edu/opendap/ORWA/contents.html                     --> Pacific North West ROMS Data
%   * http://west.rssoffice.com:8080/thredds/catalog/roms/CA300m-nowcast/catalog.html   --> Monterey Bay ROMS Data
%   * http://barataria.tamu.edu/thredds/catalog.html?dataset=txla_hindcast_agg          --> Texas-Louisiana Gulf Coast ROMS Data
%   * http://thredds.aoos.org/thredds/ghost.html                                        --> Alaska ROMS

% Information on the ROMS Data:
%   * https://www.myroms.org/wiki/Numerical_Solution_Technique

% Matlab tool box for ROMS 
%   * http://romsmatlab.tiddlyspot.com/


close all
clear all
clc

forcastDate  = datetime('today');


%% Build URL for Pacific Northwest ROMS data -----------------------------------------------------------
% * File indexing: 01-Jan-2019  --> 5114
openDap_url  = 'http://ingria.coas.oregonstate.edu:80/opendap/ORWA/ocean_his_';
forcastIndex = 5114 + days(forcastDate - datetime('01-Jan-2019'));
openDap_url  = sprintf('%s%i_%s.nc', openDap_url, forcastIndex, forcastDate);
% -----------------------------------------------------------------------------------------------------



%% URL for Monteray Bay ROMS data ------------------------------------------------------------------------------
% * This data set seems to end on Feb 23 2018??
openDap_url = 'http://west.rssoffice.com:8080/thredds/dodsC/roms/CA300m-nowcast/ca300m_das_2018022303.nc';
% -----------------------------------------------------------------------------------------------------



%% Build URL Prince William Sound --------------------------------------------------------------------------------
% forcastDate = datetime('today','Format', 'yyyy-MM-dd''T''HH:mmXXX','TimeZone','UTC');
baseDate = datetime('1970-01-01T00:00000Z','InputFormat', 'yyyy-MM-dd''T''HH:mmXXX','TimeZone','UTC');

% hr = hours( forcastDate - datetime('01-Feb-2011 06:00:00')); % This isn't working yet

lon = [166  270];                           % Acceptible Range: 0 - 449
lat = [132  193];                           % Acceptible Range: 0 - 194
depth = [0  0];                             % Acceptible Range: 0 -  15
hr = 0;

openDap_url = 'http://thredds.aoos.org/thredds/dodsC/PWS_L1_FCST.nc';
params      = sprintf('?depth[%i:1:%i],lat[%i:1:%i],lon[%i:1:%i],', depth(1), depth(2), lat(1), lat(2), lon(1), lon(2) );

temp        = sprintf('temp[%i][%i:1:%i][%i:1:%i][%i:1:%i],',   hr, depth(1), depth(2), lat(1), lat(2), lon(1), lon(2) );
salt        = sprintf('salt[%i][%i:1:%i][%i:1:%i][%i:1:%i],',   hr, depth(1), depth(2), lat(1), lat(2), lon(1), lon(2) );
u           = sprintf(   'u[%i][%i:1:%i][%i:1:%i][%i:1:%i],',   hr, depth(1), depth(2), lat(1), lat(2), lon(1), lon(2) );
v           = sprintf(   'v[%i][%i:1:%i][%i:1:%i][%i:1:%i],',   hr, depth(1), depth(2), lat(1), lat(2), lon(1), lon(2) );

zeta        = sprintf('zeta[%i][%i:1:%i][%i:1:%i],',            hr, lat(1), lat(2), lon(1), lon(2) );
time        = sprintf('time[%i]',                               hr );

openDap_url = strcat(openDap_url, params, temp, salt, u, v ,zeta, time );
%openDap_url = strcat(openDap_url, params,  time );

clear lon lat depth hr params temp salt u v zeta time 
% -----------------------------------------------------------------------------------------------------



%% Make the Call =============================================================================================================================
disp('Getting Data From OPENDaP')
disp(openDap_url)
finfo = Extract_NetCF_File(openDap_url);

clear openDap_url

if any(lon > 180)               % Longitude is in degrees East
    lon =  -(360 -lon);          % Convert into degrees West
end

% Get Geotiff
disp('Getting Geotiff')
[geoImage, geoData] = GetMapfromWMS( lon, lat );
geoshow(geoImage, geoData);

%%

[lon, lat] = meshgrid(lon,lat);

figure
hold on
geoshow(geoImage, geoData);
%scatter(lon(:),lat(:),'.','Cdata',temp(:)) 
quiver(lon ,lat ,u', v')




%% ROMS --> Show Bathymetry

time    = 1;                                        % Time index
surface = zeta(:,:,time);                           % Ocean Surface at the indexed time
param   = temp(:,:,40,time);                        % Parameter to be displayed

param(param == max(param(:)))       = NaN;          % Get rid of values in areas that are not in the ocean
surface(surface == max(surface(:))) = NaN;

Zparam = zeros(size(param))-10;

figure('Name','ROMS Stuff','NumberTitle','off')
hold on
surf(lon_rho ,  lat_rho, -h,'EdgeColor','none')                                             % Plot Bathymetry
c = contour3(lon_rho, lat_rho, -h, 40, 'k');                                                % Add contour lines to bathymetry
surf(lon_rho , lat_rho, surface,'EdgeColor','none', 'FaceColor','c', 'FaceAlpha','0.20' )  % Show the ocean surfaces 
hold off
view(0,30)


%% Display 3D Structure of a parameter

%contourLines = contourc(h,40);

contourLines = contourdata( contourc(h,40) );

dZ = max(h(:))/40;

figure
hold on
for ii = 1:40
    
    param   = temp(:,:,ii,time);                        % Parameter to be displayed
    param(param == max(param(:))) = NaN;                % Get rid of values in areas that are not in the ocean
    
    Zparam = -ones(size(param))* max(h(:)) + dZ*ii;
    
    scatter3(lon_rho(:) ,  lat_rho(:), Zparam(:) ,'.', ...
        'CData',param(:),           ...
        'MarkerEdgeAlpha', 0.20,    ...
        'MarkerFaceAlpha', 0.20);
    
end
colorbar
view(30,45)




%% Examin the surface currents

% Xvel = u(:,:,40,1);     % Size(u) = 309,  522,  40,  12
% Yvel = v(:,:,40,1);     % Size(v) = 310,  521,  40,  12

u(u == -9999.0) = NaN;
v(v == -9999.0) = NaN;

[x,y] = meshgrid(lon,lat);

figure
quiver(x,y,u(:,:,1),v(:,:,1))


%% Get the corners of the data area

max(lon_rho(:))
min(lon_rho(:))

max(lat_rho(:))
min(lat_rho(:))



% OLCI AND footprint:"Intersects(POLYGON((-122.12 49.987, -122.12 40.659, -129.99 40.659, -129.99 49.987, -122.12 49.987 )))"

