% This file imports a text file with GPS into a two column array.
% WayPoints = [Latitude, Longitude]

% Inputs:
%   filePath --> File path to the text document that is to be imported
%   fileName --> Name of the text document that is to be imported

% Outputts:
%   wayPoints --> Two column array with the imported waypoints [Latitude, Longitude] of type double

close all
clear all
clc
format long

filePath = 'C:\Users\jaanderson.FORTLEWIS\Desktop\';
fileName = 'waypoints1.txt';

    
    
    idy = 1;
    index = 1;

    file = strcat(filePath,fileName);
    
    if (exist(file))
        fprintf(' Importing waypoint list: %s\t', fileName)
    end
    
    
    logfiles = dir(fullfile(file));
    


    % Open the file, scan in the text and then close the file
    fid = fopen(strcat(filePath,logfiles(idy).name),'r');
    A = textscan(fid,'%s','Delimiter',{' '});
    fclose(fid); 

    text = A{1,1}; % Extract the charactures from their cell

    % Cycle throught each line of the text and parse the data -----------------------------------------------------------------------------------------------------------------
    [a,z] = size(text);
    for j = 1:a

        line = text{j,1};
        [z,b] = size(line);
       
        hold = {}; % Array to hold the parsed data
        char = '';
        
        try
           % Cycle through each charature in the line, seprat the values by ';' and place the values into the holding array
           hold = strsplit(line,{' ',','});
        
           % Assigh the parced values to their place in wayPoints -----------------------------------------------------------------------------------------------------
           wayPoints(index,1) = str2num(cell2mat(hold(1)));    % Latitude
           wayPoints(index,2) = str2num(cell2mat(hold(2)));    % Longitude

            index = index+1;
        catch
            
            fprintf('\n\t**Error on line %d of %s**\n',index,fileName)
           
        end % End try statement
        
    end 
   fprintf('\tDone\n\n') 
%end % End function