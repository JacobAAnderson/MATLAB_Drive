%% Batch process a directory
close all
clear all
clc


%% Choose Data Type --> Uncomment the block that you want
fileType = 'jpg';
fun = @(file) ConverPhoto(file,105);


% Load Data
ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
filePath = ui.GetDir('photos',fileType);                                      % Get folder GUI

if isempty(filePath)                                                        % Check if the Input box was cnaceled
    disp('Get Directory Input box Cancelled')
    return                                                                  % End script if there isn't an input for it to use
end

Data = BatchDir(filePath,fileType,fun);



%% Function to be applied
function conversion =  ConverPhoto(filename, imSize)

conversion = struct('Origonal', [], 'Converted', [], 'Name', '');

I  = imread(filename);

conversion.Origonal = I;

if size(I,3) > 1, I = rgb2gray(I); end

I = imresize(I,[imSize imSize]);

[path, name, ext] = fileparts(filename);

imwrite(I,fullfile(path, [name,'.jpg']));


conversion.Converted = I;
conversion.Name = [name,ext];

end