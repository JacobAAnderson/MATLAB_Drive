% This file imports a text file with GPS into a two column array.
% WayPoints = [Latitude, Longitude]

% Inputs:
%   filePath --> File path to the text document that is to be imported
%   fileName --> Name of the text document that is to be imported

% Outputts:
%   wayPoints --> Two column array with the imported waypoints [Latitude, Longitude] of type double


function wayPoints = ImportWayPoints(filePath,fileName,delimiter)
    
    idy = 1;
    index = 1;

    files = strcat(filePath,fileName);
    
    if (exist(files))
        fprintf(' Importing waypoint list: %s\t-->',files)
    
        logfiles = dir(fullfile(files));

        % Open the file, scan in the text and then close the file
        fid = fopen(strcat(filePath,logfiles(idy).name),'r');
        A = textscan(fid,'%s','Delimiter',delimiter);
        fclose(fid); 

        text = A{1,1}; % Extract the charactures from their cell

        % Cycle throught each line of the text and parse the data -----------------------------------------------------------------------------------------------------------------
        [a,z] = size(text);
        for j = 1:a

            line = text{j,1};
            %[z,b] = size(line);

            %hold = {}; % Array to hold the parsed data
            

            try
                % Cycle through each charature in the line, seprat the values by ';' and place the values into the holding array
                hold = strsplit(line,{' ',',',';','\t'});

                % Assigh the parced values to their place in wayPoints -----------------------------------------------------------------------------------------------------
                wayPoints(index,1) = str2num(cell2mat(hold(1)));    % Latitude
                wayPoints(index,2) = str2num(cell2mat(hold(2)));    % Longitude

                index = index+1;
            catch

                fprintf('\n\t**Error on line %d of %s**\n',j,fileName)

            end % End try statement

        end
    end
   fprintf('\tDone\n\n') 
end % End function