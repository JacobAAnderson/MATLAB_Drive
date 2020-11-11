
clear all
close all
clc


% Get the file-path of the directory containing the STL files
UserInputs = SavedUserInputs(mfilename);

filePath = UserInputs.getDir('STLFolder','STL');                            % Get folder GUI

if filePath == 0                                                            % Check if the Data Input box was cnaceled
    disp('Input box Cancelled.  Script Terminated')
    return                                                                  % End script if there isn't an input for it to use
end

% Read in STL files and create a single object in the desired initial orientation
files = fullfile(filePath,'*.STL*');                                        % Craete Search path for ALL .STL files in the named directory
files = dir(fullfile(files));                                               % Get names of the stl files
fprintf(' %d files found at %s \n\n',size(files,1),filePath);

if size(files,1) == 0                                                       % Exit the script if no files are found
    disp('Process Terminated')
    return
end


X = [];                                                                     % Arrays to concatinate the incoming patch objects with
Y = [];
Z = [];
C = [];

for idy = 1: size(files,1)                                                  % Read in the stl files in the directory
    fileName = fullfile(filePath,files(idy).name);
    [x, y, z, c] = stlread(fileName);
    
    X = [X,x];                                                              % Concatinate them into a single object
    Y = [Y,y];
    Z = [Z,z];
    C = [C,c];
    
end

clear x y z c idy files fileName filePath                                   % Get rid of unneeded variables


figure('NumberTitle','off','Name','Origonal Object')                        % Show what came in
axis equal
patch(X, Y, Z, C, 'FaceAlpha', 1)
xlabel('X')
ylabel('Y')
zlabel('Z')
view(45,30)



%% Tranf form arrays into the desired starting orientation

% =======================================================================================================
% Movements needed to get the object into its starting position =========================================
% Change these values!!!!!!

X_translate = 0;                                                            % Translation along X axis
Y_translate = 100.5939;                                                     % Translation along Y axis
Z_translate = 60.6;                                                         % Translation along Z axis

roll  = -90;                                                                % Rotation angle about X axis
pitch = 90;                                                                 % Rotation angle about Y axis
yaw   = 0;                                                                  % Rotation angle about Z axis
% =======================================================================================================
% =======================================================================================================


Rx = [ 1 0 0; 0 cosd(roll) -sind(roll); 0 sind(roll) cosd(roll)];           % X rotation matrix with roll  angle
Ry = [ cosd(pitch) 0 sind(pitch); 0 1 0; -sind(pitch) 0 cosd(pitch)];       % Y rotation matrix with pitch angle
Rz = [cosd(yaw) -sind(yaw) 0; sind(yaw) cosd(yaw) 0; 0 0 1];                % Z rotaiton matrix with yaw   angle

R = Rx*Ry*Rz;                                                               % Full rotation matrix

trans = [ X_translate; Y_translate; Z_translate ];                          % Translation vector

TR = [R, trans; 0 0 0 1];                                                   % Translation-Rotaion matrix

Xo = NaN(size(X));                                                          % Prealocate the size of the manipulated arrays
Yo = NaN(size(X));
Zo = NaN(size(X));


[~,objSize] = size(X);                                                      % Get the size of the graphics object

for a = 1:objSize                                                           % Perform rotaion on the Array
    
    % Extract pints from origonal patch array and apply transform
    p1 = TR * [ X(1,a); Y(1,a); Z(1,a); 1 ];
    p2 = TR * [ X(2,a); Y(2,a); Z(2,a); 1 ];
    p3 = TR * [ X(3,a); Y(3,a); Z(3,a); 1 ];
    
    % Put tranformed points in to the new patch array
    Xo(:,a) = [ p1(1); p2(1); p3(1) ];
    Yo(:,a) = [ p1(2); p2(2); p3(2) ];
    Zo(:,a) = [ p1(3); p2(3); p3(3) ];
    
end

clear a R Rx Ry Rz TR trans roll pitch yaw p1 p2 p3 X_translate Y_translate Z_translate   % Clean up unneeded variables


% Plot the new patch object
figure('NumberTitle','off','Name','Starting Configuration')
axis equal
patch(Xo, Yo, Zo, C, 'FaceAlpha', 1)                                        % Plot Object
line(xlim, [0 0], [0 0],'Color','red')                                      % Cross hairs at origin
line([0 0], ylim, [0 0],'Color','red')
line([0 0], [0 0], zlim,'Color','red')
xlabel('X')
ylabel('Y')
zlabel('Z')
view(45,30)


%% Get animation values

[fileName, filePath] = UserInputs.getFile('DataFile','txt');                % GUI to get the name and file path of a file

if fileName == 0                                                            % Check if the Data Input box was cnaceled
    disp('Input box Cancelled')
    return                                                                  % End script if there isn't an input for it to use
end

file = fullfile(filePath, fileName);
data = load(file);

time   = data(1,:);                                                         % Extract data into meanifull variables
theta  = data(2,:);
plunge = data(3,:);

clear data file fileName filePath                                                 % Clean up unneeded variables



%% Create Master Frame
clear Frames                                                                % Make sure an old instanciation of Frames in not present

% Create the frame with the figure in the starting orientation
%  * 'Position',[0 0 1920 1080] is for HD quality output
fig = figure('NumberTitle','off','Name','Animation','Position',[0 0 1920 1080]);

% Axis for the wing figure -------------------------------------------------------------------------------------------------------
ax1  = axes(fig,'units','normalized','Position',[0.0 0.2 0.7 .6]);          
patch(Xo, Yo, Zo, C, 'FaceAlpha', 1)
axis equal
ax1.View = [45,30];
drawnow

ax1.XLimMode = 'manual';                                                     % Set axes limits to manual so that the axes do not shift during the animation
ax1.YLimMode = 'manual';
ax1.ZLimMode = 'manual';

ax1.XLabel.String = 'X';
ax1.YLabel.String = 'Y';
ax1.ZLabel.String = 'Z';


% Axis for the pitch angle plot -------------------------------------------------------------------------------------------------------
ax2  = axes(fig,'units','normalized','Position',[0.67 0.6  0.3 0.3]);       % Axis for theat plot
title('Pitch Angle $\theta$ [deg]','interpreter','latex')
xlabel('Time [s]','interpreter','latex')
ylabel('$\theta$','interpreter','latex')

ax2.XLim = [min(time(:)), max(time(:))*1.1];
ax2.YLim = [min(theta(:))*1.1, max(theta(:))*1.1];

ax2.XLimMode = 'manual';                                                     % Set axes limits to manual so that the axes do not shift during the animation
ax2.YLimMode = 'manual';


% Axis for the plunge distance plot -------------------------------------------------------------------------------------------------------
ax3  = axes(fig,'units','normalized','Position',[0.67 0.1 0.3 0.3]);       % Axis for theat plot
title('Plung [m]','interpreter','latex')
xlabel('Time [s]','interpreter','latex')
ylabel('Plung distance [m]','interpreter','latex')

ax3.XLim = [min(time(:)), max(time(:))*1.1];
ax3.YLim = [min(plunge(:))*1.1, max(plunge(:))*1.1];

ax3.XLimMode = 'manual';                                                     % Set axes limits to manual so that the axes do not shift during the animation
ax3.YLimMode = 'manual';




%% Make Animation Frames

Xi = NaN(size(X));                                                          % Prealocate the sixe of the manipulated arrays
Yi = NaN(size(X));
Zi = NaN(size(X));

[~,numbData] = size(time);                                                  % Get number of animation points

for b = 1:numbData
    
    roll = 100 * theta(b);
    trans = [0;0; 1000 * plunge(b)];
    Rx =    [ 1 0 0; 0 cosd(roll) -sind(roll); 0 sind(roll) cosd(roll)];    % X rotation matrix
    TR = [Rx, trans; 0 0 0 1];
    
    for a = 1:objSize
        % Extract pints from patch array and apply transform
        p1 = TR * [ Xo(1,a); Yo(1,a); Zo(1,a); 1 ];
        p2 = TR * [ Xo(2,a); Yo(2,a); Zo(2,a); 1 ];
        p3 = TR * [ Xo(3,a); Yo(3,a); Zo(3,a); 1 ];
        
        % Assign new points to the patch array
        Xi(:,a) = [ p1(1); p2(1); p3(1) ];
        Yi(:,a) = [ p1(2); p2(2); p3(2) ];
        Zi(:,a) = [ p1(3); p2(3); p3(3) ];
        
    end
    
    % Draw new figure on frame ---------------------------------------------------------------------------------------------------------------
    cla(ax1)                                                                % Clear the previous figure
    ax1.Title.String = {'Time [s]', num2str(time(b))};                      % Update time stamp in title
    hold on
    patch(ax1,Xo, Yo, Zo, C, 'FaceAlpha', 0.09, 'EdgeAlpha', 0.05, 'FaceColor',[0.9 0.9 0.9], 'EdgeColor', [0.9 0.9 0.9]) % Ghost image of the object in its origional pose
    patch(ax1,Xi, Yi, Zi, C)                                                % New Object after transform
    hold off
    
    % Plot pitch angle -----------------------------------------------------------------------------------------------------------------------
    cla(ax2)                                                                % Clear previous data points
    plot(ax2, time(1:b), theta(1:b))                                        % Plot pitch angle
    title(ax2,  'Pitch Angle $\theta$ [deg]','interpreter','latex')         % Make sure the title and axes labels dont change
    xlabel(ax2, 'Time [s]','interpreter','latex')
    ylabel(ax2, '$\theta$','interpreter','latex')
    ax2.XLim = [min(time(:)), max(time(:))*1.1];                            % Make sure the axes limits don't change 
    ax2.YLim = [min(theta(:))*1.1, max(theta(:))*1.1];
    
    % Plot plunge height ---------------------------------------------------------------------------------------------------------------------
    cla(ax3)
    plot(ax3, time(1:b), plunge(1:b))
    title(ax3,  'Plung [m]','interpreter','latex')
    xlabel(ax3, 'Time [s]','interpreter','latex')
    ylabel(ax3, 'Plung distance [m]','interpreter','latex')
    ax3.XLim = [min(time(:)), max(time(:))*1.1];
    ax3.YLim = [min(plunge(:))*1.1, max(plunge(:))*1.1];
    
    drawnow                                                                 % Make sure image renders before the frame is captured
    
    Frames(b) = getframe(fig);                                              % Capture the figure as a movie frame
    
    
end



%% Play movie  --> Nice review step but not necesary
frameRate = 30;                                                             % Frame Rate for the movie

movFig = figure('Name','Wing Animation','NumberTitle','off','Position',[0 0 1920 1080]);  % Create figure and axes for the movie
movAxes = axes(movFig,'Units', 'normalized','Position',[0 0 1 1] );
axis off
movie(movAxes, Frames, 1, frameRate);


%% Save animation as a movie
frameRate = 30;                                                             % Frame Rate for the movie

[fileName, filePath] = UserInputs.saveFile('MoviFile','mp4');               % GUI to get file path to save file

if fileName == 0                                                            % Check if the Data Input box was cnaceled
    disp('Input box Cancelled')
    return                                                                  % End script if there isn't an input for it to use
end

videoPathName = fullfile(filePath,fileName);                                % Make full file path for the video

% Create a video from the animation
v = VideoWriter(videoPathName,'MPEG-4');                                    % Make vidoe object
v.FrameRate = frameRate;                                                    % Frame rate
open(v);
writeVideo(v,Frames);
close(v)

% Save the movie frames for later use
[framesFilepath,framesName,framesExt] = fileparts(videoPathName);
saveFrames = fullfile(framesFilepath,[framesName,'_Frames.mat']);
save(saveFrames,'Frames')
