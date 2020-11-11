% Side Scan Sonar Feature Extraction
% Jacob Anderson
% May 17, 2019


close all
clear all
clc

% Gather Data ----------------------------------------------------
ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
file = ui.getFile('image','mat');                                           % GUI to get the name and file path of a file

if isempty(file)                                                            % Check if the Data Input box was cnaceled
    return                                                                  % End script if there isn't an input for it to use
end

load(file)

if ~exist("image", 'var')                                                   % Check that the image varialbe exists in the workspace
    disp("Image was Not present in the loaded file")                        % End the script if image is not present
    return
end

image = imcrop(image);

file = ui.saveFile('image','jpeg');                                           % GUI to get file path to save file

if ~isempty(file)                                                            % Check if the Data Input box was cnaceled
    imwrite(image,file)
end

return

%% Blob Detection
close all
clc


y1 = 2*image - imdilate(image, strel('square',7));

[T, ~ , ~] = Neutrosophic_Decomposition(image);

sigma = 6;
image_filtered = imgaussfilt(T, sigma);

bw = imbinarize(image_filtered);

figure('Name', 'Binary Image')
subplot(2,1,1)
imshow(image)
subplot(2,1,2)
imshow(T) 
drawnow

hblob = vision.BlobAnalysis;

hblob.CentroidOutputPort = false;
hblob.MaximumCount = 3500;
hblob.Connectivity = 4;
hblob.MaximumBlobArea = 6500;
hblob.MinimumBlobArea = 0;
hblob.LabelMatrixOutputPort = true;
hblob.OrientationOutputPort = true;
hblob.MajorAxisLengthOutputPort = true;
hblob.MinorAxisLengthOutputPort = true;
hblob.EccentricityOutputPort = true;
hblob.ExtentOutputPort = true;
hblob.BoundingBoxOutputPort = true;

[AREA,BBOX,MAJOR,MINOR,ORIENT,ECCEN,EXTENT,LABEL] = hblob(bw);


% imshow(LABEL*2^16)

numberOfBlobs = length(AREA);
allowableAxis = (((MAJOR./MINOR) > 3.8) & (AREA > 200) & (abs(rad2deg(ORIENT))<10) & (EXTENT> 0.6));

idx = find(allowableAxis);
keeperBlobsImage = ismember(LABEL, idx);

% imshow(keeperBlobsImage)

LABEL = bwlabel(keeperBlobsImage, 4);
for i =1:length(idx)
    BBOX_OUT((i),1:4) = BBOX(idx(i),1:4);
end

NUM_BLOBS = length(idx);

figure('Name', 'Lables')
subplot(2,1,1)
imshow(LABEL*2^16)
subplot(2,1,2)
imshow(keeperBlobsImage)
drawnow





%% Neutrosophic Decomposition -------------------------------
% [T, ~ , ~] = Neutrosophic_Decomposition(image);


%% Diffusion Map
% patchDim = 10;
% configParams = [];
% 
% [patches, topleftOrigin] = im2patch(T, patchDim);
% 
% [K, nnData] = calcAffinityMat(T, configParams);
% 
% [diffusion_map, Lambda, Psi, Ms, Phi, K_rw] = calcDiffusionMap( K, configParams );
% 
% diffusionCoordIm = getDiffusionCoordIm( im, idxPatches, diffusion_mapX );
% 
% C = plotDiffusionMapinColor(diffusionCoordIm, idxPatches, diffusionMap, inds, figId, plotTitle,step);

