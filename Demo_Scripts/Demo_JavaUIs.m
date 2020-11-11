%% Demo Java Components
close all
clear all
clc

%% Standard Java JScrollBar
jScrollbar = javax.swing.JScrollBar;
jScrollbar.setOrientation(jScrollbar.HORIZONTAL);
javacomponent(jScrollbar,[10,40,200,20]);


%% Standard Java JSlider (20px high if no ticks/labels, otherwise use 45px)
jSlider1 = javax.swing.JSlider;
javacomponent(jSlider1,[10,70,200,45]);
set(jSlider1, 'Value',84, 'MajorTickSpacing',20, 'PaintLabels',true);  % with labels, no ticks


%% Vertical slider
jSlider2 = javax.swing.JSlider;
[jhSlider, hContainer] = javacomponent(jSlider2,[100,300,100,40]);
set(jSlider2, 'Value',72, ...
    'Orientation',jSlider2.VERTICAL, ...
    'MinorTickSpacing',5, ... 
    'PaintLabels',true);

set(hContainer,'position',[100,300,40,100]); %note container size change

hjSlider = handle(jSlider2, 'CallbackProperties');
hjSlider.StateChangedCallback = @(hjSlider,eventData) disp(get(hjSlider,'Value'));
% set(hjSlider, 'StateChangedCallback', @myCallback);  %alternative


%% Range Slider
jRangeSlider = com.jidesoft.swing.RangeSlider(0,255,0,255);  % min,max,low,high
jRangeSlider = javacomponent(jRangeSlider, [200,200,200,80],gcf);
set(jRangeSlider, ...
    'Opaque',false, ...
    'MajorTickSpacing',25, ...
    'MinorTickSpacing',5, ...
    'PaintTicks',false, ...
    'PaintLabels',false, ...
    'StateChangedCallback',@myCallbackFunc);




function myCallbackFunc(hjSlider,~)

fprintf('High Value: %3d \tLow Value: %3d \n\n ', hjSlider.getLowValue(), hjSlider.getHighValue() )

end