%% Import Waypoints  ==========================================================================================================================================================
function wayPoints = ImportWayPoints(filePath,fileName,delimiter)

% Add delimiter to filePath for later concatenation with file names
if ispc
    file = strcat(filePath,'\',fileName);
else
    file = strcat(filePath,'/',fileName);
end


if exist(file,'file')
    
    % Open the file, scan in the text and then close the file
    fid = fopen(file,'r');
    A = textscan(fid,'%s','Delimiter',delimiter);
    fclose(fid);
    
    text = A{1,1}; % Extract the charactures from their cell
    
    % Cycle throught each line of the text and parse the data -----------------------------------------------------------------------------------------------------------------
    index = 1;
    
    for j = 1:length(text)
        
        line = text{j,1};
        
        try
            % Cycle through each charature in the line, seprat the values by ';' and place the values into the holding array
            hold = strsplit(line,{' ',',',';','\t'});
            
            % Assigh the parced values to their place in wayPoints -----------------------------------------------------------------------------------------------------
            wayPoints(index,1) = str2num(cell2mat(hold(1)));    % Latitude
            wayPoints(index,2) = str2num(cell2mat(hold(2)));    % Longitude
            
            index = index+1;
        catch
            
            fprintf('\n\t**Error on line %d of %s**\n',j,fileName)
            
        end
        
    end
    
    
else
    disp('Waypoint file does not Exist')
    wayPoints = nan;
end

end