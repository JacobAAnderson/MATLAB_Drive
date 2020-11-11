%% Land Mask From Geotiff
%  Jacob Anderson
%  March 5 2017

% Simple color segmentation of land and water from Geotiff

function [landMask, geoImLON, geoImLAT] = LandMaskfromGeoImage( geoImage, geoData, plotResults )

[m, n, ~] = size(geoImage);
x = linspace(geoData.LongitudeLimits(1),geoData.LongitudeLimits(2),n);
y = linspace(geoData.LatitudeLimits(2), geoData.LatitudeLimits(1), m);
[geoImLON, geoImLAT] = meshgrid(x,y);


lab_he = rgb2lab(geoImage);
ab = lab_he(:,:,2:3);
ab = im2single(ab);
nColors = 3;
pixel_labels = imsegkmeans(ab,nColors,'NumAttempts',3);
BW = pixel_labels==2;
BW = ~BW;

h = fspecial('average',[10 10]);                            % Create an averaging filter to smoth edges
image = imfilter(BW,h);                                     % Applie the filter

landMask = bwareaopen(image,200);                       % Fill in the out of bounds areas
landMask = bwareaopen(~landMask,150);               % Fill in the in bounds areas
landMask = ~landMask;                               % Re-invert the ocupancy grid so that the in bounds areas are represented by 1s


if nargin > 2 && plotResults

figure('Name','Land Mask Results','Numbertitle','off')
hold on
geoshow(geoImage, geoData);
s = scatter(geoImLON(landMask), geoImLAT(landMask),'.');
hold off

title('Land Mask Results')
xlabel('Longitude')
ylabel('Latitude')
legend(s,'Land Mask')

end

end