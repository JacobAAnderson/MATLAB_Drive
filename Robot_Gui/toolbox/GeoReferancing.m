

function GeoReferancing(~)
global MainWindow;   % Main GUI window
disp('Georeferancing')

panalName = 'Geo-Ref';

MainWindow.GUIpanels(panalName,'make','position',[0.7 0 0.1 1],'backgroundcolor',[0.8510    0.3294    0.1020]);  % Create putton panel for primary GUI functions

MainWindow.Header.String = 'Select Data Input';                             % Prompt to select data

% Make buttons in the button panel --> {Buttin type, Button Name, Callback function and their inputs, Callback function for slider listeners
MainWindow.GUIbuttons( panalName,'add', ...
    {'pushbutton', 'Draw  Waypoints', @DrawWayPoints,                                                     [], '';   %  1
     'pushbutton', 'Enter Waypoints', @NumWayPoints,                                                      [], '';   %  2
     'space',      ' ',               [],                                                                 [], '';   %  3 
     'pushbutton', 'Draw  Area',      @(~,~,~) set(MainWindow.Figure,'WindowButtonDownFcn', @MouseInput), [], '';   %  4
     'pushbutton', 'Enter Area',      @NumArea,                                                           [], '';   %  5
     'space',      ' ',               [],                                                                 [], '';   %  6
     'pushbutton', 'D.R. path',       @CalcDRPath,                                                        [], '';   %  7
     'space',      ' ',               [],                                                                 [], '';   %  8
     'pushbutton', 'Done',            @GeoDone,                                                           [], '';   %  9
    });

Panel = MainWindow.Panels(panalName);                                       % Put the new version of the panel into its mapped container
Panel.button(7).Push.Visible = 'off';                                       % Hide D.R path until data is avalable

end

function DrawWayPoints(~,~)
global MainWindow DATA 
axes(MainWindow.Axes)
[DATA.Waypoints.x, DATA.Waypoints.y] = getpts(MainWindow.Axes);

MainWindow.Addlayer('Waypoints', [DATA.Waypoints.x, DATA.Waypoints.y], 'Main');                       % Add datapoints to the GUI

Panel = MainWindow.Panels('Geo-Ref');                                       % Put the new version of the panel into its mapped container
Panel.button(7).Push.Visible = 'on';                                       % Hide D.R path until data is avalable


% MainWindow.GUIbuttons( 'Geo-Ref','add', ...
%     {'pushbutton', 'D R path', {@CalcDRPath, [DATA.Waypoints.x, DATA.Waypoints.y]}, [], '';   %  1
%     });


end


function CalcDRPath(~,~)
global ui DATA

path = [DATA.Waypoints.x, DATA.Waypoints.y];

[distance, heading, time] = WayPoint2DeadReckoning(path, 2, 'knots');

time = time ./60; % Convert to minuts

filePath = ui.SaveFile('Course','txt');                                         % Get folder GUI
if isempty(filePath), return, end                                           % End script if there isn't an input for it to use

Write2File(filePath, 'Dead Reckoning Course:\n', []);
Write2File(filePath, 'Heading: %3.2f [deg], Duration: %2.3f [s], Dist: %5.2f [m]\n', [heading',time', distance']);
Write2File(filePath, '\n\nWaypoints:\n', []);
Write2File(filePath, '%f, %f\n', path);

end


function GeoDone (~,~)
global MainWindow;   % Main GUI window

MainWindow.GUIpanels('Geo-Ref','delete');  % Close the button panel for the occupancy grid tools

MainWindow.Header.String = 'Select Data Input';                             % Prompt to select data

end


function MouseInput(fig, win)
global DrawArea
persistent record


if isempty(record)
    record = false;
end

record = ~record;

if record
    DrawArea = [];
    fig.WindowButtonUpFcn = @MouseInput;
    fig.WindowButtonMotionFcn = @GetArea;
    win.Source.CurrentPoint
    
else
    fig.WindowButtonUpFcn = '';
	fig.WindowButtonDownFcn = '';
    fig.WindowButtonMotionFcn = '';
    FindColors( DrawArea )
end

end


function NumArea(~,~)

prompt = {'Enter Waypoint List: Longitude, Latitude'};
title = 'Make Area';
dims = [20 50];
definput = {['-107.846407, 37.290866';'-107.845710, 37.290768']};
answer = inputdlg(prompt,title,dims,definput);

[n,~] = size(answer{1});
if n > 2
    disp('Adding Area')
    answer
    % Do stuff

else
    disp('Not enough points to create an area')
end
    
end

function FindColors( DrawArea )
disp('find Colors')

DrawArea

end


function GetArea(~,win)
global DrawArea

DrawArea = [DrawArea; win.Source.CurrentPoint];
end

