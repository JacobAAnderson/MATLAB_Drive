%% Get Figure from movie frames
close all
clear all
clc


%% Get data to be processed
disp('Select Video')
try
    [DataFileName,DataFilePath] = uigetfile('*.mp4','Select Video',suggestions{1});
catch
    [DataFileName,DataFilePath] = uigetfile('*.mp4','Select Video');
end

if DataFileName == 0                % Check if the Data Input box was cnaceled
    disp('Input box Cancelled')
    return                          % End script if there isn't an input for it to use
end


FigureGUI(fullfile(DataFilePath, DataFileName))



%% GUI
function FigureGUI(filePath)

global axes1;
global video;

v = VideoReader(filePath);
vidWidth = v.Width;
vidHeight = v.Height;

video = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'), 'colormap',[]);

index = 1;
while hasFrame(v)
    video(index).cdata = readFrame(v);
    index = index +1;
end

nof = floor(v.Duration/v.FrameRate); %number of frames

% Make GUI to find figures
fig1 = figure('Name','Find Figure','NumberTitle','off');
axes1 = axes('Units','normalized','Position',[0,0,0.8,1]);
image(video(1).cdata,'Parent',axes1);
axis off
%movie(axes1,video,v.FrameRate)

panel = uipanel(fig1,'Title','Functions','FontSize',12,'BackgroundColor',[0.6 0.8 1], 'Units','normalized', 'Position', [0.81  0    0.2 1 ]);

slider1_txt1 = uicontrol(panel,'Style','text','String','Time',  'BackgroundColor',[0.6 0.8 1],           'Units','normalized', 'Position', [0.15 0.95 0.7  0.04]);
slider1_txt2 = uicontrol(panel,'Style','text','String','Start', 'BackgroundColor',[0.6 0.8 1],           'Units','normalized', 'Position', [0.03 0.92 0.12 0.04]);
slider1_txt3 = uicontrol(panel,'Style','text','String','End',   'BackgroundColor',[0.6 0.8 1],           'Units','normalized', 'Position', [0.77 0.92 0.12 0.04]);
slider1 = uicontrol(panel,'Style', 'slider','Min',0,'Max',v.Duration,'Value',0,'SliderStep',[1/nof 5/(nof)], 'Units','normalized', 'Position', [0.15 0.94 0.6  0.03]);
sliderListener1 = addlistener(slider1,'ContinuousValueChange',  @set_Frame);

B1 = uicontrol(panel,'Style', 'pushbutton', 'String', 'Play',      'Units','normalized', 'Position', [0.12 0.74  0.7 0.04],'Callback', @Get_Frame);
B1 = uicontrol(panel,'Style', 'pushbutton', 'String', 'Stop',      'Units','normalized', 'Position', [0.12 0.68  0.7 0.04],'Callback', @Get_Frame);
B1 = uicontrol(panel,'Style', 'pushbutton', 'String', 'Get Frame', 'Units','normalized', 'Position', [0.12 0.6   0.7 0.04],'Callback', @Get_Frame);

end

function set_Frame(source,~)
global axes1;
global video;
global frame;

frame = round(source.Value);
image(video(frame).cdata,'Parent',axes1);

end

function Get_Frame(~,~)
global video;
global frame;

figure('Name','Figure Name','NumberTitle','off');
image(video(frame).cdata);
axis off
end
