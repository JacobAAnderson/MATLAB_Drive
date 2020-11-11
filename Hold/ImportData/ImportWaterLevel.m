% This file imports a text file with GPS into a two column array.
% WayPoints = [Latitude, Longitude]



function [datum, data] = ImportWaterLevel(filePath,txtName)

% Inputs:
%   filePath --> File path to the Water level offset text document and .log files that are to be imported
%   fileName --> Name of the text document that is to be imported

% Outputts:
%   datum --> Referance point to measure water offset from
%               * Elevation in feet
%   data  --> Water level data [date, water level]
%               * Date as a datenum i.e.: days since Jan 1, 2000
%               * Water leve as an elevation in feet


offset_txt = strcat(filePath,txtName);

if (exist(offset_txt))
    fprintf(' Importing Water Offset Log: %s\t-->',offset_txt)
    
    % Open the file, scan in the text and then close the file
    fid = fopen(offset_txt,'r');
    cell = textscan(fid,'%s','Delimiter','n');
    fclose(fid);
    
    text = cell{1,1}; % Extract the charactures from their cell
    
    % Cycle throught each line of the text and parse the data -----------------------------------------------------------------------------------------------------------------
    idy = 1;
    index = 1;
    [a,z] = size(text);
    
    for j = 1:a
        
        line = text{j,1};
        
        try
            hold = strsplit(line,{' ',',',';','\t'});
            if(j == 1)
                datum = str2num(cell2mat(hold(2)));
            else
                % Assigh the parced values to their place in wayPoints -----------------------------------------------------------------------------------------------------
                date(index) = datenum(cell2mat(hold(1)),'yyyymmdd');    % Date
                waterlevel(index) = str2num(cell2mat(hold(2)));         % Recorded water levle
                index = index+1;
            end
        catch
            fprintf('\n\t**Error on line %d of %s**\n',j,offset_txt)
        end
    end
    fprintf('\tDone\n\n')
    
else
    fprintf('\n !!!! %s Does Not exist !!! \n',txtName);
    datum = nan;
    data = nan;
end

data = [date',waterlevel'];

end % End function
