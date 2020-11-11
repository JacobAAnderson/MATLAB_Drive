% Jacob Anderson
% Check Geotiff
% Oct 16, 2018

% This function checks an imported geotiff to makesure it is formated correctly.
% It will reformat the data structure if needed.

function [geoImage, geoData] = CheckGeoData(geoImage, geoData)

geoImage = im2uint8(geoImage);                                             % Convert raster image to uint8 data type for faster graphic rendering

if isempty(geoData)
    disp('Image is not a geo-referanced tiff')
    return
    
elseif ~isa(geoData,'map.rasterref.GeographicCellsReference')               % Check that geoData is in the right format for 'geoshow'
    
    geoData = georasterref('RasterSize',size(geoImage), ...                 % Reformat geoData for 'geoshow'
        'RasterInterpretation', geoData.RasterInterpretation, ...
        'LongitudeLimits',geoData.XWorldLimits, ...
        'LatitudeLimits',geoData.YWorldLimits, ...
        'ColumnsStartFrom',geoData.ColumnsStartFrom,...
        'RowsStartFrom',geoData.RowsStartFrom);
else
    disp('Geo Data checks out')
end

end