% UI Examples


%% Dialog Boxes

answer = inputdlg({'side-1', 'side-2', 'Side-3'}, 'Input data', 1, {'5','6','7'});

disp(answer)
class(answer)

s = str2double(answer);
s = sort(s);
a =s(1);
b = s(2);
c = s(3);

if (a+b) <= c
    errordlg('This Triangle Do Not Exist!!','Error!','modal')
else
    alpha = acosd((b^2+c^2-a^2)/(2*b*c));
    beta  = acosd((c^2+a^2-b^2)/(2*c*a));
    gamma = acosd((a^2+b^2-c^2)/(2*a*b));
    
    message = sprintf('The Three Angles are: %.2f, %.2f and %.2f [deg]',alpha,beta,gamma);
    msgbox(message, 'Output Data', 'modal')
end




%% UI Control with Transparent Background
close all
clear all
clc

f = figure(); % create a figure with an axes on it
ax = axes('Units','pixels', 'Position',[0 0 560 420], 'XTick',[], 'YTick',[], ...
          'Nextplot','add', 'YDir','reverse');

% read the big logo image - background of the figure
bigImage = imread('Matlab_LogoBig.png');
image(bigImage, 'parent', ax);  % Display the image in the axes

% read a smaller image - background of button
img = imread('Matlab_Logo.png');
img = imresize(img,0.25);
s = size(img);
pos = [10 10 s(2) s(1)];  %Position of the button

% Extract the portion of the image where the button will be.
F = getframe(ax,pos);  % take a screenshot of the parent figure
pb = uicontrol('Style','pushbutton', 'Units','pixels', 'Position',pos, ...
               'Callback',@(a,b)disp('push button'));
           

% as before - calculate where the button image is white.
img = double(img)/255;
index1 = img(:,:,1) == 1;
index2 = img(:,:,2) == 1;
index3 = img(:,:,3) == 1;
indexWhite = index1+index2+index3==3;

% for each pixel, replace the white portion with the parent image data
for idx = 1 : 3
   rgb = 1-img(:,:,idx);                   % To make the picture quirky change the RGB
   pImage = double(F.cdata(:,:,idx))/255;  % extract part of the image
   rgb(indexWhite) = pImage(indexWhite);   % set the white portion of the image to the parent
   img(:,:,idx) = rgb;                     % substitute the update values
end

% Update the push button image
set(pb, 'CData', img)    
    

%% UI tree
[mtree, container] = uitree('v0', 'Root','C:\', 'Parent'); % Parent is ignored
set(container, 'Parent', hPanel);  % fix the uitree Parent
