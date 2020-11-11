classdef ButtonGroup < handle
    
    properties
        Text        % Text
        Radio       % Radio buttons
        Push        % Push buttons
        Popup       % Pop up menue
        Slider      % Slider Bars
        Listener    % Active listeners for sliders
        TextEnter   % Text entry buttons
        Tag         % Tag to identify the button group
        Pos         % Verical posiiton of the button group on the GUI
        Height      % Size of the button group in y direction
        
    end
    
    methods
        
        % Constructor, Make button group -----------------------------------------------------
        function obj = ButtonGroup(parent, color, y, button)
            
            type     = button{1};       % Type of button group
            string   = button{2};       % Label for the buttons
            callBack = button{3};       % Callback function and data to pass into the call back
            listener = button{4};       % Active listener callback function for sliders
            tooltips = button{5};
            
            obj.Tag = {string};         % Tag to Identify the button group
            obj.Pos = y;                % Keep track of the group's possiotn
            
            
            % Button dimentions based on operating system
            if ispc % Dimentiond for Windows OS -----------------------------------
                  buttonHeigth = 0.05;
                    textHeight = 0.04;
                textEditHeight = 0.12;
                  sliderHeight = 0.1;
                   radioHeight = 0.04;
                
                RSHO = 0.03;            % Radio-Slider Height OffSet
                GLHO = 0.04;            % Button Group label text Heigth Offset
                TLHO = 0.01;            % Upper Text Entry Box Label Height Offset
                TEHO = 0.025;           % Upper Text Entry Box Offset
                
            else % Dimentions for Unix based systems ------------------------------
                  buttonHeigth = 0.06;
                    textHeight = 0.05;
                textEditHeight = 0.13;
                  sliderHeight = 0.11;
                   radioHeight = 0.05;
                
                RSHO = 0.04;            % Radio-Slider Height OffSet
                GLHO = 0.05;            % Slider label text Heigth Offset
                TLHO = 0.011;            % Text Entry Box Label Height Offset
                TEHO = 0.026;           % Upper Text Entry Box Offset

            end
            
            
            switch type
                
                case 'space'
                    obj.Height = buttonHeigth;                                                   % Make the hight a propertu of the button group
                    obj.Tag = 'space';

                    
                case 'text'                                                                      % Text only
                    obj.Height = textHeight;                                                     % Make the hight a propertu of the button group
                    y = y - obj.Height;                                                          % Adjust height to make room fo rthe new button
                    obj.Text = uicontrol(parent,'Style','text','units','normalized', ...         % Name the Panale as the objects parent, give object type and normalize the units
                        'Position', [0.2  y  0.6 0.04], ...
                        'BackgroundColor',color, ...
                        'HorizontalAlignment','left',...
                        'FontSize', 10, ...
                        'Tag',string, ...                                                        % Tag for callback functions to ID the source
                        'String',string );                                                       % Text to be dispalyed
                    
                
                case 'pushbutton'                                                                % Pushbuton with call back function
                    obj.Height = buttonHeigth;                                                   % Make the hight a propertu of the button group
                    y = y - obj.Height;                                                          % Adjust height to make room fo rthe new button
                    obj.Push = uicontrol(parent,'Style','pushbutton', 'units','normalized',...   % Name the Panale as the objects parent, give object type and normalize the units
                        'Position', [0.2 y 0.6 0.04], ...
                        'String', string, ...                                                    % Tag to refernece which Button the text is associated with
                        'Tag',string, ...                                                        % Tag for callback functions to ID the source
                        'Callback', callBack );                                                  % Function to execute when the buton is pushed
                    
                    if ~isempty(tooltips)
                        obj.Push.TooltipString = sprintf(tooltips);
                    end
                    
                    
                case 'popupmenu'
                    obj.Height = buttonHeigth;                                                   % Make the hight a propertu of the button group
                    y = y - obj.Height;                                                          % Adjust height to make room fo rthe new button
                    obj.Popup = uicontrol(parent,'Style','popupmenu', 'units','normalized',...   % Name the Panale as the objects parent, give object type and normalize the units
                        'Position', [0.2 y 0.6 0.04], ...
                        'String', string, ...                                                    % Tag to refernece which Button the text is associated with
                        'Callback', callBack );                                                  % Function to execute when the buton is pushed
                    
                    
                case 'edit'                                                                      % Button group for text entry box
                    obj.Height = textEditHeight;                                                 % Make the hight a property of the button group
                    y = y - obj.Height;                                                          % Adjust height to make room fo rthe new button
                    
                    % % Label for the button group __________________________________________________________________________________________
                    obj.Text{1} = uicontrol(parent,'Style','text','units','normalized', ...         
                        'Position', [0.05  y+GLHO  0.9 0.04], ...
                        'BackgroundColor',color, ...
                        'HorizontalAlignment','left',...
                        'FontSize', 9, ...
                        'Tag',string, ...                                                        % Tag for callback functions to ID the source
                        'String',string );                                                       % Text to be dispalyed
                    
                    % 1st Text Entry Box and Lable ___________________________________________________________________________________________
                    obj.Text{2} = uicontrol(parent,'Style','text','units','normalized', ...      % Lable for text entry box
                        'Position', [0.05  y+TLHO  0.2 0.04], ...
                        'BackgroundColor',color, ...
                        'HorizontalAlignment','left',...
                        'FontSize', 8, ...
                        'Tag',string, ...                                                        % Tag for callback functions to ID the source
                        'String','Max' );                                                        % Text to be dispalyed
                    
                    obj.TextEnter{1} = uicontrol(parent,'Style','edit','units','normalized', ... % Text entry box
                        'Position', [0.35  y+TEHO  0.5 0.02], ...
                        'BackgroundColor',color, ...
                        'Tag',string, ...                                                        % Tag for callback functions to ID the source
                        'Value',1, ...                                                           % Value to help identify the button 1 for max
                        'String', 'max', ...                                                     % Inital text to display
                        'Callback', callBack );                                                  % Function to execute when the buton is pushed
                    
                    % 2nd Text Entry Box and Label ___________________________________________________________________________________________
                    obj.Text{3} = uicontrol(parent,'Style','text','units','normalized', ...      % Name the Panale as the objects parent, give object type and normalize the units
                        'Position', [0.05  y  0.2 0.02], ...
                        'BackgroundColor',color, ...
                        'HorizontalAlignment','left',...
                        'FontSize', 8, ...
                        'Tag',string, ...                                                        % Tag for callback functions to ID the source
                        'String','Min' );                                                        % Text to be dispalyed
                                        
                    obj.TextEnter{2} = uicontrol(parent,'Style','edit','units','normalized', ...      % Text to Label the slider
                        'Position', [0.35  y  0.5 0.02], ...                                     % Positioned centered above the slider
                        'BackgroundColor',color, ...
                        'Tag',string, ...                                                        % Tag for callback functions to ID the source
                        'Value',0, ...                                                           % Value to help identify the button 0 for min
                        'String', 'min', ...                                                     % Tag to refernece which Button the text is associated with
                        'Callback', callBack );                                                  % Function to execute when the buton is pushed
                    
                    
                case 'slider'                                                                    % Slider button with text
                    obj.Height = sliderHeight;                                                   % Make the hight a propertu of the button group
                    y = y - obj.Height;                                                          % Adjust height to make room fo rthe new button
                    
                    % Lable for Button Group ________________________________________________________________________________________________
                    obj.Text{1} = uicontrol(parent,'Style','text','units','normalized', ...      % Text to Label the slider
                        'Position', [0.2  y+GLHO  0.6 0.04], ...                                 % Positioned centered above the slider
                        'BackgroundColor',color, ...
                        'Tag',string, ...                                                        % Tag for callback functions to ID the source
                        'String',string );
                    
                    % Label lets side of slider ______________________________________________________________________________________________
                    obj.Text{2} = uicontrol(parent,'Style','text','units','normalized', ...      % Min label
                        'Position', [0 y-0.02 0.2 0.04], ...                                     % Positioned on left side of the slider
                        'BackgroundColor',color, ...
                        'Tag',string, ...                                                        % Tag for callback functions to ID the source
                        'String','Min' );
                    
                    % Label Right Side of Slider _____________________________________________________________________________________________
                    obj.Text{3} = uicontrol(parent,'Style','text','units','normalized', ...      % Max label
                        'Position', [0.82 y-0.02 0.2 0.04], ...                                  % Positioned on right side of the slider
                        'BackgroundColor',color, ...
                        'Tag',string, ...                                                        % Tag for callback functions to ID the source
                        'String','Max' );
                    
                    % Slider and Active Listener ____________________________________________________________________________________________
                    obj.Slider = uicontrol(parent,'Style','slider','units','normalized', ...     % Slider
                        'Position', [0.2 y 0.6 0.04], ...
                        'Min',0, ...                                                             % Values are normalized
                        'Max',1, ...
                        'Value',0.5, ...                                                         % Start in the middle
                        'Tag',string, ...                                                        % Tag for callback functions to ID the source
                        'Callback', callBack);                                                   % Function to execute when the slider has been moved
                    
                    if ~isempty(listener)                                                        % Add an active listener to call a function as the slider is beeing moved
                        addlistener(obj.Slider,'ContinuousValueChange', listener);
                    end
                    
                    
                case 'sliderRadio'                                                               % Slider button with a radio pushbutton nexto it
                    obj.Height = sliderHeight;                                                   % Make the hight a propertu of the button group
                    y = y - obj.Height;                                                          % Adjust height to make room fo rthe new button
                    
                    % Lable for Button Group ________________________________________________________________________________________________
                    obj.Text{1} = uicontrol(parent,'Style','text','units','normalized', ...      % Text to Label the slider
                        'Position', [0.2  y+GLHO  0.6 0.04], ...                                 % Positioned centered above the slider
                        'BackgroundColor',color, ...
                        'String',string );
                    
                    % Label lets side of slider ______________________________________________________________________________________________
                    obj.Text{2} = uicontrol(parent,'Style','text','units','normalized', ...      % Min label
                        'Position', [0 y 0.2 0.04], ...                                          % Positioned on left side of the slider
                        'BackgroundColor',color, ...
                        'String','Min' );
                    
                    % Label Right Side of Slider _____________________________________________________________________________________________
                    obj.Text{3} = uicontrol(parent,'Style','text','units','normalized', ...      % Max label
                        'Position', [0.82 y 0.2 0.04], ...                                       % Positioned on right side of the slider
                        'BackgroundColor',color, ...
                        'String','Max' );
                    
                    % Slider and Active Listener ____________________________________________________________________________________________
                    obj.Slider = uicontrol(parent,'Style','slider','units','normalized', ...     % Slider
                        'Position', [0.2 y 0.6 0.04], ...
                        'Min',0, ...                                                             % Values are normalized
                        'Max',1, ...
                        'Value',0.5, ...                                                         % Start in the middle
                        'Tag',string);
                    if ~isempty(listener)                                                        % Add an active listener to call a function as the slider is beeing moved
                        obj.Listener = addlistener(obj.Slider,'ContinuousValueChange', listener{1});
                    end
                    
                    % Radio Button __________________________________________________________________________________________________________
                    obj.Radio = uicontrol(parent,'Style','radiobutton','units','normalized', ...  % Radio push button
                        'Position', [0.05 y+RSHO  0.1 0.04], ...                                  % Positioned left of the slider
                        'BackgroundColor',color, ...
                        'Value', 1, ...
                        'Tag',string, ...                                                        % Tag for callback functions to ID the source
                        'Callback',callBack );
                    

                    
                case 'radio'                                                                     % Radio button with text
                    obj.Height = radioHeight;                                                    % Make the hight a propertu of the button group
                    y = y - obj.Height;                                                          % Adjust height to make room fo rthe new button
                    obj.Radio = uicontrol(parent,'Style','radiobutton','units','normalized', ... % Radio push button
                        'Position', [0.2 y+0.01  0.9 0.04], ...                                  % Positioned left of the slider
                        'BackgroundColor',color, ...
                        'Value',listener , ...
                        'Tag',string, ...                                                        % Tag for callback functions to ID the source
                        'String',string,...
                        'Callback',callBack );
                    
                    
                otherwise
                    disp('Unrecognized Button Type')
                    
            end
        end
        
        
        % Disalbe the buton group ------------------------------------------------------------
        function obj = DisableButtons(obj)
            
            if ~isempty(obj.Slider) % Disable Sliders
                obj.Slider.Enable = 'off';
                
            end
            
            
            if ~isempty(obj.Push)      % Disable buttons
                obj.Push.Enable = 'off';
            end
            
            
            if ~isempty(obj.Radio)  % Disable Radio buttons
                obj.Radio.Enable = 'off';
            end
            
        end
        
        
        % Enable the button group ------------------------------------------------------------
        function obj = EnableButtons(obj)
            
            if ~isempty(obj.Slider) % Disable Sliders
                obj.Slider.Enable = 'on';
                
            end
            
            
            if ~isempty(obj.Push)      % Disable buttons
                obj.Push.Enable = 'on';
            end
            
            
            if ~isempty(obj.Radio)  % Disable Radio buttons
                obj.Radio.Enable = 'on';
            end
            
        end
        
        
        % Move the enter button group vertically to a new location ---------------------------
        function obj = MoveButtons(obj,delta)
            
            p = properties(obj);                                            % Get the properties of the class
            
            obj.Pos = delta;                                                % Get the differance in height
            
            for aa = 1: length(p)                                           % Cycle through the class propertes and delete graphics objects
                for bb = 1:length(obj.(p{aa}))
                    
                    if isa(obj.(p{aa})(bb),'matlab.ui.control.UIControl')
                        obj.(p{aa})(bb).Position = obj.(p{aa})(bb).Position + [ 0 delta 0 0];
                        
                    elseif isa((obj.(p{aa})(bb)),'cell') && isa(obj.(p{aa}){bb},'matlab.ui.control.UIControl')
                        obj.(p{aa}){bb}.Position = obj.(p{aa}){bb}.Position + [ 0 delta 0 0];
                        
                    else
                        continue
                    end
                end
            end
            
        end
        
        
        % Delete the Buttons and text insid a button group -----------------------------------
        function obj = DeleteButtons(obj)
            
            p = properties(obj);                                    % Get the properties of the class
            
            for aa = 1: length(p)                                   % Cycle through the class propertes and delete graphics objects
                
                for bb = 1:length(obj.(p{aa}))
                    
                    if isa(obj.(p{aa})(bb),'matlab.ui.control.UIControl')
                        delete(obj.(p{aa})(bb));
                        
                    elseif isa((obj.(p{aa})(bb)),'cell') && isa(obj.(p{aa}){bb},'matlab.ui.control.UIControl')
                        delete(obj.(p{aa}){bb});
                        
                    else
                        continue
                    end
                end
            end
        end
        
    end
end