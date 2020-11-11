% Difustion map


close all
clear all
clc

ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class

% Get file with file path GUI ----------------------------------------------------------------------------------------------
file = ui.getFile('image','mat');                                            % GUI to get the name and file path of a file

if isempty(file)                                                            % Check if the Data Input box was cnaceled
    return                                                                  % End script if there isn't an input for it to use
end

load( file )

% Kernal Function
kernal = @(x,y,s) exp( -(x-y).*(x-y)./ s );

sigma = 0.5;

d = zeros(size(image));

for ii = 1:numel(image)
    
   d_sub = kernal(ii, image, sigma);
    
   d(ii) = sum( d_sub(:) );
       
end