% This file imports a text file with GPS into a two column array.
% WayPoints = [Latitude, Longitude]

% Inputs:
%   filePath --> File path to the Water level offset text document and .log files that are to be imported
%   fileName --> Name of the text document that is to be imported

% Outputts:
%   bathy --> Four column array with the GPS data, adjusted water depth in meters, and the calculated offset in meters [lon, lat, depth, offset]


function [bathy, DATA] = WaterLevelCorrection(filePath,txtName)

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
                date(index) = datenum(cell2mat(hold(1)),'yyyymmdd');    % Latitude
                waterlevel(index) = str2num(cell2mat(hold(2)));    % Longitude
            end
            index = index+1;
        catch
            fprintf('\n\t**Error on line %d of %s**\n',j,offset_txt)
        end
    end
    fprintf('\tDone\n\n')
    
else
    fprintf('\n !!!! %s Does Not exist !!! \n',txtName);
    bathy = nan;
end

% DATA.depth(1) --> Total water Column [m]
% DATA.date     --> Date: mm.dd
DATA = LogConcat(filePath);

fprintf('\n Applying Offsets \n\n')

offsetDates = datenum({DATA.date},'mm.dd,yyyy');

rawDepth = cell2mat({DATA.depths}');   % [total water column, vehicle depth] in meters

GPS = fliplr(cell2mat({DATA.GPS}'));

% __Syntax__ y = fixpt_interp1(xdata,ydata,x,xdt,xscale,ydt,yscale,rndmeth)
elevation = fixpt_interp1(date,waterlevel,offsetDates,'double',2^-14,'double',2^-14,'Nearest');

offset = (datum - elevation).*0.3048;   % Calculate offset in meters

Depth = rawDepth(:,1) + offset;

bathy = [GPS,Depth,rawDepth,offset];

end % End function
