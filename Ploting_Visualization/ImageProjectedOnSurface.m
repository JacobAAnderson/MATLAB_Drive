% Jacob Anderson
% Project image onto a 3D surface
% October 15, 2018

close all   % Fresh start
clear all
clc


% Check for existing input preferances
UserInputs = SavedUserInputs(mfilename);

% Get Geo-tiff -----------------------------------------------------
disp('Select Geo-tif')
try
    [geotiff_name, geotiff_filepath] = uigetfile('*.tif*','Select Geo-tiff',UserInputs.Preferences('Geotiff'));
catch
    [geotiff_name, geotiff_filepath] = uigetfile('*.tif*','Select Geo-tiff');
end

if geotiff_name == 0                                                        % Check if the Data Input box was cnaceled
    disp('No Getiff selected')                                              % End script if there isn't an input for it to use
else
    geotiff_fullfilepath = fullfile(geotiff_filepath, geotiff_name);        % Get full file path
    UserInputs.Preferences('Geotiff') = geotiff_fullfilepath;               % Save file path
    UserInputs.SavePreferences();
    
    [geoImage, geoData] = geotiffread(geotiff_fullfilepath);                % Read Goetiff
    
    geoData = CheckGeoData(geoImage, geoData);


end



imshow(geoImage)

%% Create Arrays for the surface
[m,n,~] = size(geoImage);
x = -n/2:n/2;
y = -m/2:m/2;

[X,Y] = meshgrid(x,y);

Z = cos(10*X/n).*cos(5*Y/m);

%% Plot Surface with image

I = rot90(geoImage);   % Rotate the image 180 degrees to display with the botom of the image forward
I = rot90(I);

figure('Name', 'Surface with image projected on top of it');
surf(X, Y, Z, I,'EdgeColor','none'); % Plot surface


