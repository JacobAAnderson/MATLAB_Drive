%% Ecomapper data to Graylog

% Jacob Anderson
% Robotic GNOME Laboratory
% Fort Lewis College
% Durango, CO, 81301

% 2017

% All Rights Reserved.

% Contact: Ryan N. Smith - rnsmith@fortlewis.edu

% This is open software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version. It would be appreciated to provide
% credit to the originating source when appropriate.
%
%  This is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You can find the GNU General Public License at <http://www.gnu.org/licenses/>.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc


% Find path to this file
[thisFile_path,thisFile_name,~] = fileparts(mfilename('fullpath'));

% Create an file for saving data or acces previously saved data
UserInputs = SavedUserInputs(thisFile_name);


%% Get File path to log files
disp('Select Data folder');                  % Prompt to select data
try
    filePath = uigetdir(UserInputs.Preferences('DataFolder'),'Select Data Folder');
catch
    filePath = uigetdir('*.tif*','Select Geo-tiff');
end

if filePath == 0                                        % Check if the Data Input box was cnaceled
    disp('Input box Cancelled')
    return                                              % End script if there isn't an input for it to use
elseif  ~exist(filePath,'dir')
    error(' --- Invaled file path to data ---')
else
    fprintf('\n Importing EcoMapper logs... \n\n')
    UserInputs.Preferences('DataFolder') = filePath;    % Add file path to saved inputs
    UserInputs.SavePreferences() 
end


%% Cycle through all the .log files ---------------------------------------------------------------------------------------------------------------
files = fullfile(filePath,'*IVER2*.log');
logfiles = dir(fullfile(files));
fprintf(' %d files found at %s \n\n',size(logfiles,1),files);


for idy = 1: size(logfiles,1)
    tic
    fprintf(' * Checking file %d of %d: %s \n',idy,size(logfiles,1),logfiles(idy).name)
    
    % Read Log file ==============================================================================
    % Check that the files have data. i.e.. more than just a header
    if (logfiles(idy).bytes < 950 )
        fprintf('\t--> Empty file, skip \n')
        
        % If the file has data, read it in line by line and parse the data
    else
        fprintf('\t--> Has Data, pulling in \n')
        
        try
            fileID = fopen(fullfile(filePath,logfiles(idy).name),'r');
            dataArray = textscan(fileID,'%s', 'HeaderLines', 1 ); %, 'Delimiter', delimiter, 'ReturnOnError', false);
            fclose(fileID);
        catch exception
            
            warning off backtrace
            warning('Error Reading log file \n\t\t%s \n\t\t%s \n\t\tError in: %s, Line: %d \n\n\t\tSkipping over the file \n\n',exception.identifier,exception.message,exception.stack(1).name,exception.stack(1).line)
            warning on backtrace
            continue
        end
    end
    
    dataArray = dataArray{1};
    
    %% Send data to Graylog =============================================================================
    
    
 %   header = jsonencode(containers.Map( ...
 %       {'Content-Type', 'Accept'}, ...
 %       {'application/json', 'application/json'} ))
    
    body = jsonencode(containers.Map( ...
        { 'version','host','short_message','full_message', 'timeStamp'}, ...
        { '1.1', 'Particle.IO', 'Data from Sensor Buoy 1', '{{PARTICLE_EVENT_VALUE}}','{{PARTICLE_PUBLISHED_AT}}'})) 

    
    url = 'http://logger.fortlewis.edu:12202/gelf';     
    
    options = weboptions('RequestMethod', 'post');
    options.MediaType = 'application/json';
    options.Timeout = 60;
    options.HeaderFields = {'Content-Type', 'application/json'; ...
                            'Accept', 'application/json' }
    
    response = webwrite(url, body, options)
    
    
  %  for line = 1: size(dataArray,1)    
  %  end
    
end



