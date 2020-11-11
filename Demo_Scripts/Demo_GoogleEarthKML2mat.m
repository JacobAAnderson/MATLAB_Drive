
close all
clear all
clc


ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(true);                                                     % Indicate whether new data should be selescted by the user


% Get file with file path GUI ----------------------------------------------------------------------------------------------
file = ui.GetFile('Data','kml');                                            % Get file path GUI
if isempty(file), return, end                                               % End script if there isn't an input for it to use

% Data = GoogleEarthKML2mat(file);


%% Export as csv

filePath = ui.SaveFile('Course','csv');                                     % Get folder GUI
if isempty(filePath), return, end                                           % End script if there isn't an input for it to use

[dir, ~, ~] = fileparts(filePath);

for data = GoogleEarthKML2mat(file)
    
    name = data.name;
    
    name = regexprep(name, ' ', '_');
    
    path = [data.coordinates(:,1), data.coordinates(:,2)];
    
    [distance, heading, time] = WayPoint2DeadReckoning(path, 0.8, 'm/s');
    
    time = time ./60; % Convert to minuts
    
    filePath = fullfile(dir,[name,'.txt']);
   
    Write2File(filePath, 'Mission:%s\n\n', name);
    Write2File(filePath, 'Dead Reckoning Course:\n', []);
    Write2File(filePath, 'Heading: %-3.2f [deg], Duration: %-2.3f [min], Dist: %-5.2f [m]\n', [heading',time', distance']);
    Write2File(filePath, '\n\nWaypoints:\n', []);
    Write2File(filePath, '%f, %f\n', path);
    
    %     writematrix(data.coordinates, ['/Users/jake/Desktop/',name,'.txt'])
    
end



