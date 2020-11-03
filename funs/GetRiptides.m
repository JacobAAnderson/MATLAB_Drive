% Get Riptides
% Jacob Anderson
% RDML, OSU, Corvallis OR.
% Sept 22, 2020



% Get Riptide -------------------------------------------------------------
function [RT, geotiff, bathy] = GetRiptides(ui, num, names, filters)

% Get file Paths
path(num).data        = '';                                                 % File Paths
path(num).acomms      = '';
path(num).pAcommsHand = '';
path(num).waypoints   = '';
path(num).offset      = '';

for ii = 1:num
    
    path(ii).data = ui.GetDir( sprintf('%s Data', names{ii} ), 'csv');           % Get file path GUI
    if isempty(path(ii).data), return, end                                       % End script if there isn't an input for it to use
    
    path(ii).acomms = ui.GetDir( sprintf('%s Acomms', names{ii} ), 'csv');       % Get file path GUI
    if isempty(path(ii).acomms), return, end
    
    path(ii).pAcommsHand = ui.GetDir(sprintf('%s pAcomms', names{ii} ),'txt');   % GUI to get file path to pAcommsHandler logs
    if isempty(path(ii).pAcommsHand), return, end 

    path(ii).offset = ui.GetFile(sprintf('%s Off Set', names{ii} ), 'txt');      % GUI to get Get water level offset
    
    path(ii).waypoints = ui.GetDir( sprintf('%s Waypoints', names{ii} ), 'txt');  % Get file path GUI
    if isempty(path(ii).waypoints), return, end
    
end

geotiff = ui.GetGeoTiff;                                                    % Get Geotiff
if isempty(geotiff.Image), return, end


bathFile = ui.GetFile('bathymetry','mat');                                  % Get file path GUI
if isempty(bathFile), return, end

bathy = LoadBathymetry(bathFile);                                           % Load Bathymetry Map

% Instanciate Riptide Data

for ii = num:-1:1
    
    if nargin >= 4
        RT(1:ii) = Riptide_Data(path(ii).data, path(ii).offset, names{ii}, filters);         % Create Object with Riptide Data data logs and give vehicle name
        RT(ii) = RT(ii).Add_Acomms(path(ii).acomms, filters);               % Add acoustic communications logs
    else
        RT(1:ii) = Riptide_Data(path(ii).data, names{ii});                  % Create Object with Riptide Data data logs and give vehicle name
        RT(ii) = RT(ii).Add_Acomms(path(ii).acomms);                        % Add acoustic communications logs
    end
    
    
    RT(ii) = RT(ii).Add_pAcommsHandler(path(ii).pAcommsHand);               % Add pAcommsHandler logs
    RT(ii) = RT(ii).Add_Waypoints(path(ii).waypoints);                      % Add mission waypoints
%     RT(ii) = RT(ii).Get_Manifest("skipidel","alt & acomms");                % Create a table summerising the data and acomms logs, skip instances where the vehicle is idel and disply the table

end


for rt = RT
    rt.Disp_Manifest;                                                       % Print manifest
end

end

