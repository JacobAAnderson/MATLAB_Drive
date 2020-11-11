%% Start Using Data from Copernicus

% OLCI AND footprint:"Intersects(POLYGON((-121.774 37.045, -121.774 36.40, -122.749 36.40, -122.749 37.045, -121.774 37.045 )))"


clear all
close all
clc

% Get File Name
UserInputs = SavedUserInputs(mfilename);                                    % Instantiate the Saved User Inputs Class
% Get file with file path GUI ----------------------------------------------------------------------------------------------
file = UserInputs.getFile('Data','nc');                                    % GUI to get the name and file path of a file

if isempty(file)                                                            % Check if the Data Input box was cnaceled
    disp('Get File Input box Cancelled')
    return                                                                  % End script if there isn't an input for it to use
else
    fprintf('\nFile Path and name: %s \n\n', file );
end

finfo = Extract_NetCF_File(file);

clear file



%%

romsArea = [-122.12, 49.987; -122.12, 40.659; -129.99, 40.659; -129.99, 49.987; -122.12, 49.987];

lat = double(latitude)  .* double( finfo.Variables(2).Attributes(3).Value );
lon = double(longitude) .* double( finfo.Variables(3).Attributes(3).Value );

in = inpolygon(lon, lat, romsArea(:,1), romsArea(:,2));

lat = lat(in);
lon = lon(in);
alt = altitude(in);

%%
scatter3(lon, lat, alt,'.','Cdata',alt)

