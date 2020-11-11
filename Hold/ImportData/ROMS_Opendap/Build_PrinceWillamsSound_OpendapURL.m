% Build URL quire string for Prince Willams Sound Alaka
% Jacob Anderson
% March 9, 2019


% This script pulls data in from a ROMS model for the Prince Williams Sound in Alaska.
% It also retrives a geo referanced image of the area.
% Tune the parameters in the data sellection area to select the data that you want.

% NOTE: ROMS variables will be created in the workspace.

% ROMS Resources ===============================================================================================
% Info Page:
%   * http://thredds.aoos.org/thredds/ghost.html

% Information on the ROMS Data:
%   * https://www.myroms.org/wiki/Numerical_Solution_Technique


% Acceptible Data Ranges________________________________________________________________________________________
% Time in 3hr incraments: 0 - 17947  -->  01-Feb-2011 06:00:00 to 29-Aug-2017 06:00:00 
% Longitude:              0 -   449  -->  -154.5000 to -139.5333  [deg.dec]
% Latitude:               0 -   194  -->    55.0000 to   61.4667  [deg.dec]     [155    179]
% Depth:                  0 -    15  -->     0.0000 to  300.0000  [m]
%_______________________________________________________________________________________________________________

function [openDap_url, refDate] = Build_PrinceWillamsSound_OpendapURL(hrs, lon, lat, depth)


if any(hrs > 17947), hrs(hrs > 17947) = 17947; end
if any(hrs < 0),     hrs(hrs < 0)     =     0; end

if any(lon > 449), lon(lon > 449) = 449; end
if any(lon < 0),   lon(lon < 0)     = 0; end

if any(lat > 194), lat(lat > 194) = 194; end
if any(lat < 0),   lat(lat < 0)     = 0; end

if any(depth > 15), depth(depth > 15) = 15; end
if any(depth < 0),  depth(depth < 0)  =  0; end

% % Maybe get an upto date forcast??? -------------------
% forcastDate = datetime('today');
% hr = hours( forcastDate - datetime('01-Feb-2011 06:00:00') ); % This isn't working yet
% ???? ------------------------------------------------


% Referance Date for the data
refDate = datetime('1970-01-01T00:00000Z','InputFormat', 'yyyy-MM-dd''T''HH:mmXXX','TimeZone','UTC');

openDap_url = 'http://thredds.aoos.org/thredds/dodsC/PWS_L1_FCST.nc';
params      = sprintf('?depth[%i:1:%i],lat[%i:1:%i],lon[%i:1:%i],', depth(1), depth(2), lat(1), lat(2), lon(1), lon(2) );

temp        = sprintf('temp[%i:1:%i][%i:1:%i][%i:1:%i][%i:1:%i],',  hrs(1), hrs(2), depth(1), depth(2), lat(1), lat(2), lon(1), lon(2) );
salt        = sprintf('salt[%i:1:%i][%i:1:%i][%i:1:%i][%i:1:%i],',  hrs(1), hrs(2), depth(1), depth(2), lat(1), lat(2), lon(1), lon(2) );
u           = sprintf(   'u[%i:1:%i][%i:1:%i][%i:1:%i][%i:1:%i],',  hrs(1), hrs(2), depth(1), depth(2), lat(1), lat(2), lon(1), lon(2) );
v           = sprintf(   'v[%i:1:%i][%i:1:%i][%i:1:%i][%i:1:%i],',  hrs(1), hrs(2), depth(1), depth(2), lat(1), lat(2), lon(1), lon(2) );

zeta        = sprintf('zeta[%i:1:%i][%i:1:%i][%i:1:%i],',           hrs(1), hrs(2), lat(1), lat(2), lon(1), lon(2) );
time        = sprintf('time[%i:1:%i]',                              hrs(1), hrs(2) );

openDap_url = strcat(openDap_url, params, temp, salt, u, v ,zeta, time );

end