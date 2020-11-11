%% Demo support vector machine Process

close all
clear all
clc

%% Get Trainning Data
truthmap = imread('/home/jake/MATLAB-Drive/HomeWork/SDM_Project/blob1.gif');

mapSize = numel(truthmap);

numPoints = ceil(mapSize/10);

samplePoints = randi(mapSize, [1,numPoints] );

[Xt,Yt] = ind2sub(size(truthmap), samplePoints(20:end) );

trainingData = truthmap( samplePoints(20:end) );

figure
subplot(1,2,1)
ax = gca();
scatter(ax, Xt,Yt,'.');
%freezeColors(ax);
hold(ax, 'on');
imh = imshow(logical(truthmap));
hold(ax, 'off')
uistack(imh, 'bottom')

subplot(1,2,2)
scatter(Xt,Yt,'.','cdata',uint8(trainingData))



%% Train Model
Mdl = fitrsvm([Xt',Yt'],double(trainingData), 'KernelFunction', 'gaussian');


%% Make Preditions
[Xp,Yp] = ind2sub(size(truthmap), samplePoints(1:19) );

truthData = truthmap( samplePoints(1:19) );

predictedData = predict(Mdl,[Xp',Yp']);

figure
subplot(2,2,1)
ax = gca();
scatter(ax, Xp,Yp,'.');
%freezeColors(ax);
hold(ax, 'on');
imh = imshow(logical(truthmap));
hold(ax, 'off')
uistack(imh, 'bottom')

subplot(2,2,2)
scatter(Xp,Yp,'.','cdata',uint8(predictedData))

subplot(2,2,4)
scatter(Xp,Yp,'.','cdata',uint8(truthData))



