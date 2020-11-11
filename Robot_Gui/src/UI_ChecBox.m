% UI Check Box Imput


function indx = UI_ChecBox(name, list)
global indx

list_size = numel(list);
box_height = min([1/(list_size + 2), 0.06]);
box_width  = 0.75;

indx = true(size(list));


fig = figure('name',sprintf("% Selection Box",name),'numbertitle','off');

cbx(list_size) =  uicontrol(fig,'style','checkbox', ...
                                'units', 'normalized', ...
                                'Position',[0.1  (0.95 - box_height)  box_width  box_height]);

for ii = 1: list_size
    
    cbx(ii) =  uicontrol(fig, 'style','checkbox', ...
                              'String', list{ii}, ...
                              'Value', 1,...
                              'units', 'normalized', ...
                              'Position',[0.1  (0.95 - ii * box_height)  box_width  box_height], ...
                              'Callback', @CheckBox_Callback, ...
                              'UserData', ii );
end

uicontrol(fig, 'style','pushbutton', ...
               'String', 'OK',...
               'units', 'normalized', ...
               'Position',[0.4, 0.1, 0.2, box_height],...
               'Callback', {@OKButton_Callback, fig} );

drawnow;

waitfor(fig)

end


function CheckBox_Callback( ui, ~)
global indx

indx(ui.UserData) = ui.Value; 

% assignin('caller','indx',ind)

end


function OKButton_Callback(~,~,h), close(h); end
