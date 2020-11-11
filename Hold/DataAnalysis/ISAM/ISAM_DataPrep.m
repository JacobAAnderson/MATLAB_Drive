% Side Scan Sonar SLAM Data Prep
% Jake Anderson
% 9/26/2019

close all
clear all
clc


ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(true);                                                      % Indicate whether new data should be selescted by the user

% Get Ecomapper log, SSS logdoc, bathymetry map and Geotiff
emFile = ui.GetFile('EcomapperData','log');                                 % Get Ecomapper .log file
if isempty(emFile), return, end                                             % Check if the Data Input box was cnaceled

sssFilePath = ui.GetDir('SideScanSonarData','logdoc');                      % Get Side Scan Sonar .logdoc file
if isempty(sssFilePath), return, end                                        % Check if the Input box was cnaceled

bathFile = ui.GetFile('BathymetryMap','mat');                               % Get Bathymetry Map .mat file
if isempty(bathFile), return, end                                           % Check if the Data Input box was cnaceled

geotiff = ui.GetGeoTiff;                                                    % Gui to get Geotiff

% Load Data
fprintf("\n\n--- Loading Data ---\n\n")
emData     = EM_Data(emFile);                                               % Load ecomapper data
sssData    = SSS_Data(sssFilePath);                                         % Load side Scan Sonar Data
bathymetry = load(bathFile);                                                % load Bathymetry map
bathymetry = bathymetry.bathy_map;                                          % Get rid of loading structue

clear emFile sssFilePath bathFile                                           % Clean up the workspace

% Data Smoothing
fprintf("\n\n--- Filtering Data ---\n\n")
station_ID = 9410079;                                                       % NOAA station ID for Catlina IS
emData = emData.TideCorrection(station_ID);                                 % Get the Tide Level that correspondes to the when the data was colected 
emData = emData.FilterAltimiter(2);                                         % Low Pass Filter for Raw Altimiter Data, Cut data points that excced 2 std
emData = emData.FilterHeading(1);                                           % Low Pass Filter for Raw Compass Data, Cut data points that excced 1 std
emData = emData.GetVehiclePaths;                                            % Generate the vehicle's path in lat-lon, utm and dead reckoning
emData = emData.Altimeter_Path;                                             % Generate the Bathymetry readings in lat-lon, utm, dead reckoning

sssData = sssData.MakeSoanrRanges(30);                                      % Establish the ranges of the sonar returns using a max range of 30
sssData = sssData.Index2EM_Data(emData.RawData.timeStamp);                  % Index the Side Scan sonar data to the ecomapper data
sssData = sssData.GeoMesh_From_Bathymetry(bathymetry);                      % Use bathymetry to create a geomesh for the sonar track

% Prepar Data for U-Net Preditions
fprintf("\n\n--- Prepar Data For U-Net CNN ---\n\n")

savefile = ui.SaveFile('ImagePatch_for_CNN','jpeg');                        % Where to save image pathces for CNN feature extractor
if isempty(savefile), return, end 

windowSize = 512;                                                           % Size of the image patch that will be sent to the CNN --> This will be downsampled by 1/2 when saved as a jpeg
overlap    = 0.75;                                                          % Percentage of window overlap
sssData = sssData.DiceSonarTrack(windowSize, overlap, savefile);            % Dice the Entire sonar tracke into image pathces for the CNN

clear windowSize overlap savefile station_ID                                % Clean up the workspace


%% -- Get U-Net Predictions --
filePath = ui.GetDir('CNN_Prediction_Binaries','jpeg');                     % Where to get the predition images produced by the U-net CNN
if isempty(filePath), return, end 

sssData = sssData.GetDicedBinaries(filePath);                               % Get Binary Image preditions from the CNN

imageSize = 210;                                                            % Size of the image patches --> this will be down sampled by 1/2 when the image is saved
maxSize = imageSize^2 * 0.75;                                               % Maximum feature size
minSize = imageSize^2 * 0.1;                                                % Minimum feature size
sssData = sssData.GetFeaturesFromBinaries(maxSize, minSize, imageSize);     % Process Binar Images for features

% Send features to Siamese Network for matching
savefile = ui.SaveFile('ImagePatch_for_SiameseNetwork','jpeg');             % Where to Save the image pathces for the siamese network
if isempty(savefile), return, end 

sssData = sssData.ExportFeatureImages(savefile);                            % Save features' image patches as jpegs

clear imageSize maxSize minSize filePath savefile                           % Clean up the workspace


%% Get Siames Predictions
file = ui.GetFile('Siamese_FeatureMatchingData','txt');                     % Where to get Siamese network preditions for feature matching
if isempty(file), return, end 

sssData = sssData.GetFeatureMatchingData(file);                             % Get matching data and create probability matrix

savefile = ui.SaveFile('ISAM_Data','mat');                                  % Where to Save the processed data
if isempty(savefile), return, end 

save(savefile, 'emData', 'sssData', 'geotiff', 'bathymetry', '-v7.3')       % Save all the data

clear file savefile                                                         % Clean up the workspace

fprintf("\n\nDone!!!!\n\n")
