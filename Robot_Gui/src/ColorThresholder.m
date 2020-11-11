function ColorThresholder

global geoImage;

disp('Color Threshold App')

screensize = get( groot, 'Screensize' );
screenWidth  = screensize(3);
screenHeigth = screensize(4);

boxColor = [0.9020, 0.9020, 0.9020];

fig = figure('Name','Color Thresholder', ...
    'units','pixel',...
    'Position',[screenHeigth*0.1  screenWidth*0.1   screenWidth*0.8  screenHeigth*0.7],...
    'Color',[1 1 1],...
    'NumberTitle','off',...
    'MenuBar','none',...
    'ToolBar','none', ...
    'Resize','off');


ImAax = axes(fig,'Units','pixel', ...
    'Position',[fig.Position(3)*0.05 fig.Position(4)*0.15 fig.Position(3)*0.5 fig.Position(4)*0.7], ...
    'Box','off', ...
    'XColor', fig.Color, ...
    'YColor', fig.Color, ...
    'XTick',[], ...
    'YTick',[], ...
    'ButtonDownFcn',@GetArea);

image(geoImage, 'parent', ImAax);  % Display the image in the axes


% Make color select bodes
ax(3) = gobjects;
y = 0.4;
color = ['r','g','b'];

for ii = 1:3
    
    ax(ii) = axes(fig,'Units','pixel', ...
        'Position',[fig.Position(3)*0.65 fig.Position(4)*y fig.Position(3)*0.3 fig.Position(4)*0.15], ...
        'Nextplot','add', ...
        'XTick',[], ...
        'YTick',[], ...
        'XColor', boxColor, ...
        'YColor', boxColor, ...
        'Color', boxColor, ...
        'PlotBoxAspectRatio',[1 0.25 1]);
    
    
    histogram(ax(ii),geoImage(:,:,ii),255, ...
        'Normalization','countdensity',...
        'EdgeColor','non', ...
        'FaceColor',color(ii));
    
    
    jRangeSlider = com.jidesoft.swing.RangeSlider(0,255,0,255);  % min,max,low,high
    jRangeSlider = javacomponent(jRangeSlider, [fig.Position(3)*0.65 fig.Position(4)*(y-0.01) fig.Position(3)*0.3 fig.Position(4)*0.01], fig);
    set(jRangeSlider, ...
        'Opaque',false, ...
        'MajorTickSpacing',25, ...
        'MinorTickSpacing',5, ...
        'PaintTicks',false, ...
        'PaintLabels',false, ...
        'StateChangedCallback',{ @ColorSliderCallBack, color(ii), ImAax});
    
    y = y + 0.2;
    
end


% uicontrol(fig,'Style','pushbutton', ...
%     'units','pixel',...
%     'Position', [fig.Position(3)*0.75 fig.Position(4)*0.3 fig.Position(3)*0.06 fig.Position(4)*0.05], ...
%     'String', 'Choos Area', ...
%     'Callback', @(~,~) set(fig,'WindowButtonDownFcn', @MouseInput) );

% t = timer;
% t.ExecutionMode = 'fixedRate';                  % Exicute TimerFcn at a fixed rate
% t.TimerFcn = @PlotBW;                    % Function to be exicuted by the timer
% t.Period = 0.5;                                   % Execution Period in seconds
% t.TasksToExecute = 10;                          % Number of exicutions before stopping
% t.StopFcn = @(~,~) disp('Timer has Stopped');   % Tell the user that the timer has stopped

uicontrol(fig,'Style','pushbutton', ...
    'units','pixel',...
    'Position', [fig.Position(3)*0.75 fig.Position(4)*0.2 fig.Position(3)*0.06 fig.Position(4)*0.05], ...
    'String', 'Done', ...
    'Callback',{@LoadBW, fig} );


end

function ColorSliderCallBack(hjSlider,~, channel, ImAax)
persistent redIn greenIn blueIn
global geoImage;
global BW;


% Assign Values to RGB limits if this is the first time the function has been called
if isempty(redIn)
    [row, column, ~] = size(geoImage);
    redIn = true(row, column);
    greenIn = true(row, column);
    blueIn = true(row, column);
    
end

% Up Date color thresholds
switch channel
    case 'r'
        redIn   = (geoImage(:,:,1) <= hjSlider.getHighValue()) & ...
                  (geoImage(:,:,1) >= hjSlider.getLowValue());
        
    case 'g'
        greenIn = (geoImage(:,:,2) <= hjSlider.getHighValue()) & ...
                  (geoImage(:,:,2) >= hjSlider.getLowValue());
        
    case 'b'
        blueIn  = (geoImage(:,:,3) <= hjSlider.getHighValue()) & ...
                  (geoImage(:,:,3) >= hjSlider.getLowValue());
        
    otherwise
        disp('ColorThresholder /  ColorSliderCallBack : unrecognized RGB input')
        
end

BW = redIn & greenIn & blueIn;
PlotBW(ImAax);
end

function PlotBW(ImAax)
persistent rowGrid columnGrid
global geoImage;
global BW;


if isempty(rowGrid)
    [row, column, ~] = size(geoImage);
    ROW = 1:row;
    COLUMN = 1:column;
    [rowGrid, columnGrid] = meshgrid(COLUMN,ROW);
    
end

out = ~BW;
axes(ImAax)
imshow(geoImage)
hold on
s = scatter(rowGrid(out),columnGrid(out),out(out),'k');
s.MarkerFaceColor = 'none';
s.MarkerEdgeAlpha = 0.6;
hold off
drawnow limitrate

end

function LoadBW(~,~,fig)
disp('Color Thresholder is Done')
% stop(t);
% delete(t);
delete(fig);

end

% ========================================================================================================
% 
% function MouseInput(fig, win)
% global DrawArea
% persistent record
% 
% 
% if isempty(record)
%     record = false;
% end
% 
% record = ~record;
% 
% if record
%     DrawArea = [];
%     fig.WindowButtonUpFcn = @MouseInput;
%     fig.WindowButtonMotionFcn = @GetArea;
%     win.Source.CurrentPoint
%     
% else
%     fig.WindowButtonUpFcn = '';
% 	fig.WindowButtonDownFcn = '';
%     fig.WindowButtonMotionFcn = '';
%     FindColors( DrawArea )
% end
% 
% end
% 
% 
% function FindColors( DrawArea )
% disp('find Colors')
% 
% DrawArea
% 
% end
% 
% 
% function GetArea(~,win)
% global DrawArea
% 
% DrawArea = [DrawArea; win.Source.CurrentPoint];
% end
% 
