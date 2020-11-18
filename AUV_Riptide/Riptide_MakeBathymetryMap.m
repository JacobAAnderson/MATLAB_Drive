% Demo Bathymetry_Map
% Jacob Anderson
% 9/16/2019

close all
clear all
clc


% Get Data Sources
ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(false);                                                     % Indicate whether new data should be selescted by the user

geotiff = ui.GetGeoTiff;                                                    % Gui to get Geotiff ---> (lon, lat) is optional for locating file on your computer but needed for retriving map from Google / web map server



%% Bathymetry Map from Riptide and Lutra

% --- Get data from Riptides ---
paths.rt = ui.GetDir('RT_Data','csv');                                      % Get folder GUI
if isempty(paths.rt), return, end                                           % End script if there isn't an input for it to use

paths.offset = ui.GetFile('Off Set','txt');  

paths.Lutra = ui.GetDir('LU_Data','csv');                                   % Get folder GUI
if isempty(paths.Lutra), return, end      


% --- Get Data From Riptides ---
filters = {'GPS_fix', 0};                                                    % Data filter {"Data_Field", value to be filterred out}
%           'ALT_ALTITUDE', 0};  

RT = Riptide_Data(paths.rt, paths.offset, '21x', filters);                  % Read in Riptide Files
rt_bathy = cat(1, RT.rawData.bathymetry);


% [alt, idx] = FilterRawAltimeter(rt_bathy(:,3), 16);                         % Filter Altimiter Data
% rt_bathy = [rt_bathy(idx,1:2), -alt(idx)];


rt_bathy(:,3) = -rt_bathy(:,3);                                             % Invert depth measurments


% --- Get Data From Lutra ---- 
lu_data = BatchDir(paths.Lutra,'csv',@Get_LU_Data);                         % Read in Lutra data Files

data = cat(1, lu_data.data);                                                % Combine data sources
data = unique(data,'rows');                                                 % Get Rid of duplicate data points

 
% -- Data is expected to be [lat, lon, depth] --
data = [data(:,3), data(:,2), data(:,4)];

data = [data; rt_bathy];                                                    % Combine riptide and Lutra data


% --- clean Up the data ----
out = any(isnan(data),2) | ...
      any(isinf(data),2) | ...
      data(:,3) == 0     | ...
      data(:,3) > 16;
      

data(out,:) = [];

% perimiter = csvread(file4);
% in = inpolygon(data(:,1), data(:,2), perimiter(:,2), perimiter(:,1)); 
% data = data(in,:);


clear file1 file2 filePath tmp data1 data2 data3 rt_data rt_bathy lu_data



%% Make the Bathymetry Map
bathy_map = Bathymetry_Map;                                                 % Instanciate bathymetry map
bathy_map = bathy_map.MakeMap(data);                                        % Make Bathymetry map from rando data

bathy_map.PlotMap(geotiff);

% bathy_map = bathy_map.MapSmoothing;

bathy_map.PlotMap(geotiff);

% Make Bathymetry map with GP Regresstion --> still being developed
% bathy_gp = Bathymetry_Map;
% bathy_gp = bathy_gp.MakeMap_GP(data, perimiter);




%% Display Bathymetry
% fig1 = bathy_map.Plot_3DModel(250, 30);                                     % Plot Bathymetry Model
% 
% if ~isempty(geotiff.Image)
%     fig2 = bathy_map.PlotMap(geotiff);                                      % Plot Map over the Geotiff 
% end


%% Save Map
ui = ui.NewData(true);                                                      % Indicate that new data should be selescted by the user
saveFile = ui.SaveFile('Bathymetry','mat');                                 % GUI to get file path to save variable / object
if isempty(saveFile), return, end                                           % End script if there isn't an input for it to use
disp("Saving Bathymetry Map")
save(saveFile, 'bathy_map')



%% Other functions
function Data = Get_LU_Data(file)
fprintf(' * Importing Riptide Data Log: %s\n',file)
Data.data = csvread(file);

end



function [alt, idx] = FilterRawAltimeter(alt, threshold)

idx = false(size(alt));

if nargin == 2
    idx(alt > threshold) = true;
end

idx(alt >= 0)   = true;
idx(isnan(alt)) = true;
idx(isinf(alt)) = true;

alt(idx) = 0;

if nargout == 2
    idx = ~idx;
end

A = 1;              % Process Model
B = 1;              % Control Model
C = 1;              % Observation Matrix

Q = 0.1;             % Uncertainty in the process model
R =std(alt);        % Measurment Standard Deviation

x = 0;              % State
sig = R;            % Error in the state
u = 0;              % Control input

for ii = 2:numel(alt)
    
    if alt(ii) == 0
        sig = R * 2;
        continue 
    end
    
    z = alt(ii);
    
    [x, sig] = Kalman_Filter(x, sig, u, z, A, B, C, Q, R);
    
    alt(ii) = x;
end


end






