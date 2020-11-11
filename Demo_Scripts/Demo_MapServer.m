% Demo Map server


%% Mapp Server
close all
clear all
clc

lonLimit = [-120.0 , -90.0];
latLimit = [   0.0,   90.0];

nasa = wmsfind('nasa', 'SearchField', 'serverurl');
layer = nasa.refine('bluemarbleng',  'SearchField', 'layername', 'MatchType', 'exact');
[geoImage, geoData] = wmsread(layer,'Latlim',latLimit,'Lonlim',lonLimit);
geoshow(geoImage, geoData)

%% Hi-res Imagry
close all
clear all
clc

numberOfAttempts = 5;
attempt = 0;
info = [];

%serverURL = 'http://raster.nationalmap.gov/arcgis/services/Orthoimagery/USGS_EROS_Ortho_NAIP/ImageServer/WMSServer?request=GetCapabilities&service=WMS';
 serverURL = 'http://raster.nationalmap.gov/arcgis/services/Orthoimagery/USGS_EROS_Ortho_1Foot/ImageServer/WMSServer?';
while(isempty(info))
    try
        info = wmsinfo(serverURL,'TimeoutInSeconds',180);
        orthoLayer = info.Layer(1);
    catch e 
        
        attempt = attempt + 1;
        if attempt > numberOfAttempts
            throw(e);
        else
            fprintf('Attempting to connect to server:\n"%s"\n', serverURL)
        end        
    end
end
