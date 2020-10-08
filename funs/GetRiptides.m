% Get Riptides
% Jacob Anderson
% RDML, OSU, Corvallis OR.
% Sept 22, 2020



% Get Riptide -------------------------------------------------------------
function [RT, geotiff, bathy] = GetRiptides(ui, num, names, filters)

% Get file Paths
path_data{num}        = '';                                                 % File Paths
path_acomms{num}      = '';
path_pAcommsHand{num} = '';
path_waypoints{num}   = '';

for ii = 1:num
    
    path_data{ii} = ui.GetDir( sprintf('Data%d',ii), 'csv');                % Get file path GUI
    if isempty(path_data{ii}), return, end                                  % End script if there isn't an input for it to use
    
    path_acomms{ii} = ui.GetDir( sprintf('Acomms%d',ii), 'csv');            % Get file path GUI
    if isempty(path_acomms{ii}), return, end
    
    path_pAcommsHand{ii} = ui.GetDir(sprintf('pAcomms%d',ii),'txt');                      % GUI to get file path to pAcommsHandler logs
    if isempty(path_pAcommsHand{ii}), return, end 

    path_waypoints{ii} = ui.GetDir( sprintf('waypoints%d', ii), 'txt');    % Get file path GUI
    if isempty(path_waypoints{ii}), return, end
    
end

geotiff = ui.GetGeoTiff;                                                    % Get Geotiff
if isempty(geotiff.Image), return, end


bathFile = ui.GetFile('bathymetry','mat');                                  % Get file path GUI
if isempty(bathFile), return, end

bathy = LoadBathymetry(bathFile);                                           % Load Bathymetry Map

% Instanciate Riptide Data

for ii = num:-1:1
    
    if nargin >= 4
        RT(1:ii) = Riptide_Data(path_data{ii}, names{ii}, filters);         % Create Object with Riptide Data data logs and give vehicle name
        RT(ii) = RT(ii).Add_Acomms(path_acomms{ii}, filters);               % Add acoustic communications logs
    else
        RT(1:ii) = Riptide_Data(path_data{ii}, names{ii});                      % Create Object with Riptide Data data logs and give vehicle name
        RT(ii) = RT(ii).Add_Acomms(path_acomms{ii});                            % Add acoustic communications logs
    end
    
    
    RT(ii) = RT(ii).Add_pAcommsHandler(path_pAcommsHand{ii});               % Add pAcommsHandler logs
    RT(ii) = RT(ii).Add_Waypoints(path_waypoints{ii});                      % Add mission waypoints
%     RT(ii) = RT(ii).Get_Manifest("skipidel","alt & acomms");                % Create a table summerising the data and acomms logs, skip instances where the vehicle is idel and disply the table

end


for rt = RT
    rt.Disp_Manifest;                                                       % Print manifest
end

end

