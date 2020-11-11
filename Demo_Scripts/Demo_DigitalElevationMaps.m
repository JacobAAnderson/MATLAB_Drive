
% Digital Bathymetry data

% NOAA Digital Elivation Maps --> Small areas high resolution  
%   * https://www.ngdc.noaa.gov/mgg/coastal/coastal.html

% General Bathymettic Chart of the Ocean --> Large area low resolution
%   * https://www.gebco.net/

% Prince William Sound, Alaska 8 Arc-second MHHW Coastal Digital Elevation Model
%   * https://data.noaa.gov//metaview/page?xml=NOAA/NESDIS/NGDC/MGG/DEM/iso/xml/638.xml&view=getDataView&header=none



%% Get Alaska 8 Arc-second MHHW Coastal Digital Elevation Model 
openDap_url ='https://www.ngdc.noaa.gov/thredds/dodsC/regional/prince_william_sound_8_mhhw_2009.nc';



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


%% ShoW Alaska Coastal Digital Elevation Model

Band1(Band1 == -99999 ) = NaN;


[X,Y] = meshgrid(lon,lat);

scatter3(X(:), Y(:), Band1(:),'.', 'Cdata', Band1(:) );
