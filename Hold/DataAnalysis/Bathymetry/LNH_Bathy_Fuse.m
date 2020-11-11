% Bathymetry Map from Mat
% Jacob Anderson
% 3/04/2020

close all
clear all
clc


%% Choose Data 

fileType = 'mat';
fun = @(file) load(file);


% Load Data
ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(false);                                                     % Indicate whether new data should be selescted by the user

filePath = ui.GetDir('Data',fileType);                                      % Get folder GUI

if isempty(filePath), return, end                                                        % Check if the Input box was cnaceled

Data = BatchDir(filePath,fileType,fun);


if isempty(Data)                                                            % Check that data is present
    disp("Empty Data Structure")
    disp("Doulbe check the directroy that you chose")
    return                                                                  % End script if there isn't any data to process
end


for ii = 1:size(Data,1)
    
    bathy = Data(ii).Bathymetry;
    lat = Data(ii).LAT;
    lon = Data(ii).LON;
    
    lat(isnan(bathy)) = [];
    lon(isnan(bathy)) = [];
    bathy(isnan(bathy)) = [];
    
    Data(ii).Bathymetry = bathy(:);
    Data(ii).LON = lon(:);
    Data(ii).LAT = lat(:);
    

end

data(:,3) = cat(1, Data.Bathymetry);
data(:,2) = cat(1, Data.LON);
data(:,1) = cat(1, Data.LAT);

clear bathy lon lat filePath Data



%% Make Bathymetry Map
bathy = Bathymetry_Map;                                                 % Instanciate Bathymetry Map
bathy = bathy.MakeMap(data);



%% Display Bathymetry

bathy.Elevation = -bathy.Elevation;

fig = bathy.Plot_3DModel(250, 30);                                     % Plot Bathymetry Model



%% Save Map
saveFile = ui.SaveFile('Bathymetry','mat');                                 % GUI to get file path to save variable / object
if isempty(saveFile), return, end                                           % End script if there isn't an input for it to use

save(saveFile, 'bathy')


