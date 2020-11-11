%% Extract all info from a netcd file

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

varnames = Extract_NetCF_File(file);
