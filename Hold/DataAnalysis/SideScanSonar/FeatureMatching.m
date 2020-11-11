% Feature Matching
% Jacob Anderson
% May 7, 2019

close all
clear all
clc

ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class


%% Get Data
file1 = ui.getFile('Data1','mat');                                           % GUI to get the name and file path of a file

if isempty(file1)                                                           % Check if the Data Input box was cnaceled
    return                                                                  % End script if there isn't an input for it to use
end


file2 = ui.getFile('Data2','mat');                                           % GUI to get the name and file path of a file

if isempty(file2)                                                           % Check if the Data Input box was cnaceled
    return                                                                  % End script if there isn't an input for it to use
end

image1 = load(file1);
image1 = image1.image;

image2 = load(file2);
image2 = image2.image;


clear file1 file2

%% Find Mathching Features
close all
clc

% s = size(image2,2) - size(image1,2);

%image2 = image2(:, s:end);

% Prefilter image
% sigma = 6;
% image1_filtered = imgaussfilt(image1, sigma);
% image2_filtered = imgaussfilt(image2, sigma);


% Fined Features ----------------------------------------------------------
% points1 = detectHarrisFeatures(image1_filtered);
% points2 = detectHarrisFeatures(image2_filtered);

points1 = detectKAZEFeatures(image1, 'Diffusion','region', 'NumOctaves', 4); %, 'NumScaleLevels', 4, 'Threshold', 0.0000001);  % Find Blobles
points2 = detectKAZEFeatures(image2, 'Diffusion','region', 'NumOctaves', 4); %, 'NumScaleLevels', 4, 'Threshold', 0.0000001); 

% points1 = detectSURFFeatures(image1_filtered);
% points2 = detectSURFFeatures(image2_filtered);

% points1 = detectBRISKFeatures(image1_filtered);
% points2 = detectBRISKFeatures(image2_filtered);

% points1 = detectFASTFeatures(image1_filtered); %,'MinContrast',0.3);
% points2 = detectFASTFeatures(image2_filtered); %,'MinContrast',0.3);

% points1 = detectMSERFeatures(image1_filtered);
% points2 = detectMSERFeatures(image2_filtered);


% Eliminate features in sonar traugh --------------------------------------
midLine = floor(size(image1,1)/2);
ii = 1;
while ii <= points1.Count
    
    if midLine + 50 < points1.Location(ii,1) && points1.Location(ii,1) > midLine - 50
        points1(ii,:) = [];
%         points1.Metric(ii,:) = [];
%         points1.Scale(ii,:) = [];
%         points1.Orientation(ii,:) = [];
    end
    
    ii = ii+1;
end


% Eliminate weak features -------------------------------------------------
% points1 = selectStrongest(points1,20);
% points2 = selectStrongest(points2,20);


% Extract Features --------------------------------------------------------
[features1,valid_points1] = extractFeatures(image1,points1);
[features2,valid_points2] = extractFeatures(image2,points2);


% Match Features ----------------------------------------------------------
indexPairs = matchFeatures(features1,features2, 'Method', 'Approximate', 'Unique', true); %, 'MatchThreshold', 5.75, 'MaxRatio', 1, 'Metric', 'SAD');

matchedPoints1 = valid_points1(indexPairs(:,1),:);
matchedPoints2 = valid_points2(indexPairs(:,2),:);

% Dispaly features --------------------------------------------------------
figure('Name', "Features",'NumberTitle','off')
subplot(2,1,1)
imshow( image1 )
hold on
plot(valid_points1)
hold off

subplot(2,1,2)
imshow( image2 )
hold on
plot(valid_points2)
hold off


figure; 
showMatchedFeatures( image1, image2, matchedPoints1, matchedPoints2 );
