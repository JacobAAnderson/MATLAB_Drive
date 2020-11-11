classdef (Sealed = true) GUIwindow < handle
    
    properties
        Figure            % GUI window
        Axes              % Plotting axes i.e. dispaly area
        Header            % Text Box on top of the GUI window
        Panels            % Tool bars
        Layers            % Plotting layers
        mouseCounter      % Count the turns of the mouse scrole wheel
    end
    
    properties (Access = private)
        LayerTransparancy % Variable to hold the transparancy value of the layers
        AlphaMax
        LayerOnOff        % Variable to indicate whether a layer sould be ploted
        Xgrid             % Variable to hold the Longitude mesh
        Ygrid             % Variable to hold the Latitude mesh
        geoData           % Variable to Geographic referance data from geotiff
        geoImage          % Variable to hold raster image from geotiff
        bkgnImage         % Back ground Image
    end
    
    
    methods
        
        % Make GUI window with ploting axes and button panel ---------------------------------------------------------------------------
        function obj = GUIwindow(varargin)
            % varargin: 1 cell, GUI window name --> Optional
            
            
            figName = 'GUI Window';         % Default window name in case it is omitted
            
            if nargin == 1                  % Get the name of the GUI Windos
                figName = varargin{1};
            elseif nargin > 2
                warning('Too many Inputs')
            end
            
            % Create the GUI window with a plotting area and a header
            obj.Figure = figure('Name',figName, 'units','normalized',...
                'Position',[0.1  0.1   0.7  0.7],...
                'Color',[0.8  0.8  0.8],...
                'NumberTitle','off',...
                'MenuBar','none',...
                'ToolBar','none'); %, ...
            % 'CloseRequestFcn', @obj.CloseRequest);
            
            c = uicontextmenu;
            uimenu(c,'Label','Save Figure','Callback',@obj.MakeFigure);

            
            obj.Axes   = axes(obj.Figure, ...
                'Position',[0.1  0.1   0.65 0.7], ...
                'Color',[0.25, 0.25, 0.25], ...
                'UIContextMenu', c);
            
            obj.Header = uicontrol(obj.Figure,'Style','text','units','normalized', ...
                'Position',[0  0.95  0.8  0.05], ...
                'BackgroundColor',[0.3020  0.7490  0.9294], ...
                'String', 'Hello', ...
                'FontName', 'FixedWidth', ...
                'FontSize', 14, ...
                'FontWeight','bold', ...
                'HorizontalAlignment', 'left' );
            
            
            % Setup the rest of the class properties
            obj.Panels            = containers.Map('KeyType', 'char' ,'ValueType','any');     % Tool bars
            obj.Layers            = containers.Map('KeyType', 'char' ,'ValueType','any');     % Plotting layers
            obj.LayerTransparancy = containers.Map('KeyType', 'char' ,'ValueType','double');  % Variable to hold the transparancy value of the layers
            obj.LayerOnOff        = containers.Map('KeyType', 'char' ,'ValueType','logical'); % Variable to indicate whether a layer sould be ploted
            obj.Xgrid             = [];                                                       % Variable to hold the Longitude mesh
            obj.Ygrid             = [];                                                       % Variable to hold the Latitude mesh
            obj.geoData           = [];                                                       % Variable to Geographic referance data from geotiff
            obj.geoImage          = [];                                                       % Variable to hold raster image from geotiff
            obj.mouseCounter      = 1;                                                        % Variable to indicate the number of mouse scrole wheel turns
            obj.AlphaMax          = 0.25;                                                     % Maximum Layer Trancparancy
        end
        
        
        % Manage Button Panels ------------------------------------------------------------------------------------------------------------
        function obj = GUIpanels(obj, panelName, direction, varargin)
            % varargin: variable input, background color, pose, etc
            
            switch direction
                
                case 'make'                                     % Create a new tool bar
                    
                    pose =nan;
                    backGroundColor = [0.6510    0.6510    0.6510];
                    
                    if nargin > 3               % Get panel position if available, otherwis use a default posisiton
                        
                        for ll = 1:length(varargin)
                            
                            if isa(varargin{ll},'char')
                                
                                switch lower(cell2mat(varargin(ll)))
                                    
                                    case 'backgroundcolor'
                                        backGroundColor = varargin{ll+1};
                                        
                                    case 'position'
                                        pose = varargin{ll+1};
                                        
                                    otherwise
                                        continue
                                end
                                
                            else
                                continue
                                
                            end
                        end
                    end
                    
                    
                    
                    allPanels = keys(obj.Panels);               % Get the names of all of the tool bars
                   
                    if contains(allPanels, panelName)           % If the panel in questions already exists, then exit the function
                        return
                        
                    elseif length(allPanels) >= 3                % If three panels are already present, close the third "Tool Bar" Panel
                        
                        for aa = 1: length(allPanels)
                            Panel = obj.Panels(allPanels{aa});
                            
                            x = Panel.panel.Position;
                            
                            if x(1) <= 0.79
                                obj.GUIpanels(allPanels{aa},'delete');
                            end
                            
                        end
                        
                    elseif length(allPanels) > 1
                        obj.Axes.Position = [0.1 0.1 0.55 0.7]; % Resize the GUI's dispaly area to accomidate for the new panel1
                        
                    end
                    
                    
                    if isnan(pose)
                        if isempty(allPanels)
                            pose = [0.8 0 0.2 1];
                            
                        elseif length(allPanels) == 1
                            pose = [0.8 0 0.2 0.49];
                            
                        else
                            pose = [0.7 0 0.1 1];
                            
                        end
                    end
                    
                    
                    Panel.panel = uipanel(obj.Figure,'Title', panelName,'Position',pose,'BackgroundColor',backGroundColor,'BorderType', 'etchedout','FontSize',14);
                    obj.Panels(panelName) = Panel;
                    
                    
                case 'empty'                            % Delete the contents of a panel without deleting the panel itself
                    
                    Panel = obj.Panels(panelName);      % Get the graphics objects ( buttons and text ) from the mapped container
                    names = fieldnames(Panel);          % Generate a list of the fields present in the structure
                    in = ~ismember(names,'panel');      % Remove the panel from the list of objects to delete
                    names = names(in);
                    
                    for ii = 1:length(names)            % Cycle through the graphics objects and delete them
                        % Panel.(names{ii});
                        delete(Panel.(names{ii}))
                    end
                    
                    % Reassign the generic structure to the panel's mapped container
                    obj.Panels(panelName) = struct('panel',Panel.panel,'button',gobjects,'Slider',gobjects,'Text',gobjects,'Radio',gobjects);
                    
                    
                case 'delete'                           % Delete the contents of a panel and panel itself
                    Panel = obj.Panels(panelName);      % Get the graphics objects ( buttons and text ) from the mapped container
                    delete(Panel.panel);                % Delet Graphics objec
                    remove(obj.Panels,panelName);       % Remove the panel from the mapped container
            end
            
        end
        
        
        % Add buttons to GUI window panel -------------------------------------------------------------------------------------------------
        function obj = GUIbuttons(obj, panalName, dirction, varargin)
            % varargin: array of varying size:
            % If "adding" buttons varargin is a single cell containing a N x 4 cell array of:
            % {'Button type', 'Label String', 'Callback function and values to pass to callbak function', 'Active Listener Callbak function' }
            
            % If "removing" buttons varargin is a single cell containing the name of th ebutton group to remove
            % Buttons underneath the removed button group will automaticaly be moved up
            
            Panel = obj.Panels(panalName);          % Extract the button panel from its mapped container to access the data structure
            Color = Panel.panel.BackgroundColor;    % Get panel's background color to assign to the buttons and text.
            
            switch dirction                                         % Setup the indexing variables for creating new buttons or adding more buttons
                
                case 'add'                                          % Adding buttons, determin how many objects are present and start the indes one above that.
                    y = 0.98;                                       % Starting Vertical possition of the new objects
                    if isfield(Panel, 'button')                     % Get number of buttons present on the button panel
                        numButtons = length(Panel.button) +1;
                        
                        for jj = 1:length(Panel.button)
                            y = y - Panel.button(jj).Height;
                        end
                        
                    else
                        numButtons = 1;
                        
                    end
                    
                    % Creat buttons from the top of the panel downward -----------------------------------------------------------------------------------
                    value = varargin{:};
                    [row, ~] = size(value);
                    for ii = 1: row
                        
                        Panel.button(numButtons) = ButtonGroup(Panel.panel,Color,y,value(ii,:));
                        y = y - Panel.button(numButtons).Height;                                    % Decreas the vertical location of the "button"
                        numButtons = numButtons + 1;
                        
                    end
                    
                    
                case 'remove'
                    
                    % Remove the named button
                    groups = cat(1,Panel.button.Tag);           % Get names of the button groups present
                    index = ismember(groups, varargin{1});      % Find its place in the structure
                    
                    delta = Panel.button(index).Height;         % Get size of the space being created by removing the button
                    
                    Panel.button(index).DeleteButtons();        % Delete the buttons
                    Panel.button(index) = [];                   % Remove class from from the structure
                    
                    % Move buttons up
                    pos = find(index);                          % Get Location of deleted button as an integer
                    
                    for ii = pos: length(Panel.button)          % Cycle through the buttons under the deleted button and move them up
                        Panel.button(ii).MoveButtons(delta);
                    end
                    
                    
                otherwise
                    disp('Unreckognized Button command')
            end
            
            
            obj.Panels(panalName) = Panel;                      % Return new values to the mapped container
        end
        
        
        % Set Axes Background on GUI axes -------------------------------------------------------------------------------------------------
        function obj = BackGround(obj,varargin)
            % varargin: array of variable size.
            % * Could contain: Longitude and Latatude grids or geotiff variables
            
            for i=1:length(varargin)
                try
                    switch lower(cell2mat(varargin(i)))
                        
                        case 'lon'                               % Get Window size
                            obj.Xgrid = varargin{i+1};
                            
                        case 'lat'                               % Get Window size
                            obj.Ygrid = varargin{i+1};
                            
                        case 'geodata'                               % Get Window size
                            obj.geoData = varargin{i+1};
                            
                        case 'geoimage'                               % Get Window size
                            obj.geoImage = varargin{i+1};
                            
                        case 'image'
                            obj.bkgnImage = varargin{i+1};
                        
                        otherwise
                            continue
                    end
                    
                catch
                    continue
                    
                end
            end
            
            if ~isempty( obj.geoData)
                
                obj.Axes.XLimMode = 'manual';               % Set Axis Limit managment to manual
                obj.Axes.YLimMode = 'manual';
                
                obj.Axes.XLim = obj.geoData.LongitudeLimits;    % Set new axes limits
                obj.Axes.YLim = obj.geoData.LatitudeLimits;
                
            end
            
            xlabel('Easting');                          % Lable Axes
            ylabel('Northing');
            
            obj.Plot;                                   % Plot
            
        end
        
        
        % Plot layers on GUI axes ----------------------------------------------------------------------------------------------------------------
        function obj = Plot(obj)
            
            allLayers = keys(obj.Layers);               % Get a list of all the layers present
            
            axes(obj.Axes)                              % Set the focuse on the GUI's display window
            
            
            % Draw Background image
            if isempty(obj.geoImage) && isempty(obj.bkgnImage)
                cla                                     % Clear the display area for new plots
            elseif ~isempty(obj.geoImage)
                geoshow(obj.geoImage, obj.geoData);     % Dispaly geo-tiff
            else
                imshow(obj.bkgnImage)
            end
            
            hold on                                     % Keep the new plots from being cleared as layers are added
            
            % Draw Data Layers
            for ii = length(allLayers):-1:1             % Cycle through the layers and plot them
                
                if obj.LayerOnOff(allLayers{ii})        % Asses if the layer has been turned on or off by the user
                    
                    switch allLayers{ii}
                        
                        case 'OcupancyGrid'                                 % Plot ocupancy grid
                            plotPoints = obj.Layers('OcupancyGrid');        % Get the plotting points from the mapped container
                            plotPoints = logical(plotPoints);
                            
                            if ~any(any(plotPoints)) || numel(obj.Xgrid) ~= numel(plotPoints)  % Check that there are points present
                                continue
                            end
                               
                            lon = obj.Xgrid(plotPoints);                    % Get coresponding geographic locations "geo-refrenceing arrays"
                            lat = obj.Ygrid(plotPoints);
                            plotPoints = plotPoints(plotPoints);            % Reshape plotting ploint array to match the geo-refrenceing arrays
                            
                            s = scatter(lon, lat, plotPoints,[ 0, 0.4510, 0.7412]);     % Plot as a point cloud, color --> light blue
                            s.MarkerFaceColor = 'none';                                 % Turn off the Face Color for better rendering
                            s.MarkerEdgeAlpha = obj.LayerTransparancy(allLayers{ii});   % Adjust transparancy to users settings
                            
                            
                        case 'Shape'                                        % Plot Shape file
                            plotPoints = obj.Layers('Shape');
                            plot(plotPoints(:,1),plotPoints(:,2),'b');      % Plot waypoints as a line
                            
 
                        case 'Boundaries'                                   % Plot boundarys as a line
                            plotPoints = obj.Layers('Boundaries');
                            plot(plotPoints(:,1),plotPoints(:,2),'y');
                            
                        
                        case 'Regions'                                      % Plot regions as point clouds with different colors
                            % * This function is more for design and debugging purposes.
                            Regions = obj.Layers('Regions');
                            
                            if max(Regions(:)) < 1  || numel(obj.Xgrid) ~= numel(Regions)
                                continue
                            end
                            
                            for jj = 1: max(Regions(:))
                                
                                if any(any(Regions == jj))
                                    plotPoints = ones(size(Regions));
                                    plotPoints = plotPoints(Regions == jj);
                                    lon = obj.Xgrid( Regions == jj);
                                    lat = obj.Ygrid( Regions == jj);
                                    
                                    s = scatter(lon, lat, plotPoints);
                                    s.MarkerFaceColor = 'none';
                                    s.MarkerEdgeAlpha = obj.LayerTransparancy(allLayers{ii});
                                end
                                
                            end
                            
                            
                        case 'DataPoints'                                   % Plot the goegraphic position of the loaded data points
                            
                            dataPoints = obj.Layers('DataPoints');
                            s = scatter(dataPoints(:,2), dataPoints(:,1),'g','.');
                            s.MarkerFaceColor = 'none';
                            s.MarkerEdgeAlpha = obj.LayerTransparancy('DataPoints');
 
                            
                        case 'Path'                                         % Plot the goegraphic position of the loaded data points
                            
                            dataPoints = obj.Layers('Path');
                            s = scatter(dataPoints(:,2), dataPoints(:,1),'g','.');
                            s.MarkerFaceColor = 'none';
                            s.MarkerEdgeAlpha = obj.LayerTransparancy('Path');
                            plot(dataPoints(1,2), dataPoints(1,1),'d','color','r')


                        case 'Marker'
                            dataPoints = obj.Layers('Marker');
                            if isempty(dataPoints)
                                continue
                            end
                            
                            plot(dataPoints(:,1), dataPoints(:,2),'+m')
                        
                            
                            
                        case 'Waypoints'
                            dataPoints = obj.Layers('Waypoints');
                            if isempty(dataPoints)
                                continue
                            end
                            
                            plot(dataPoints(:,1), dataPoints(:,2),'+m')
                            plot(dataPoints(:,1), dataPoints(:,2),'b')
                            plot(dataPoints(1,1), dataPoints(1,2),'om')
                            
                            
                        otherwise
                            dataPoints = obj.Layers(allLayers{ii});
                            [c,h] = contourf(obj.Xgrid, obj.Ygrid, dataPoints, 'k');
                            clabel(c,h,'FontSize',8)
                            caxis(obj.Axes,[min(dataPoints(:)), max(dataPoints(:))])
                            colorbar;
                            
                    end
                    
                else
                    continue
                end
            end
            
            drawnow limitrate       % Make plots render right away an limit the the number of updates
            hold off
            
        end
        
        
        % Manage Layers
        % Add a ploting layer with transparancy buttons -----------------------------------------------------------------------------------
        function obj = Addlayer(obj,layerName, layerData, fromPanel )
            
            if ~any(ismember(keys(obj.Panels),'Layers'))                        % Make a button panel for controling layers if it doesn't exist
                
                Panel = obj.Panels(fromPanel);                                % Resize the main panel which will be in the way
                Panel.panel.Position = [0.8 0.51 0.2 0.5];
                
                obj = obj.GUIpanels('Layers','make','position',[0.8 0, 0.2, 0.5],'backgroundcolor',[0.4940, 0.1840, 0.5560]); % Make the "Layers" Panel
                %Panel = obj.Panels('Layers');
                %Panel.panel.BackgroundColor = [0.6 0.7 1];
            end
            
            obj.Layers(layerName) = layerData;                                  % Add Layer to the mapped container
            
            
            % Check to see if this is a new layer, as opposed to a layer that is being recreated.
            allLayers = keys(obj.Layers);
            allIO     = keys(obj.LayerOnOff);
            
            newLayer = ~ismember(allLayers,allIO);
            if any(newLayer)                                                    % If it is a new layer, create the coresponing transparance and On-off buttons and varialbes
                obj.LayerOnOff(layerName) = true;
                obj.LayerTransparancy(layerName) = 0.25;
                
                % Add button for that layer
                % {'Button type', 'Label String', 'Callback function and values to pass to callbak function', 'Active Listener Callbak function' }
                obj.GUIbuttons( 'Layers','add', ...
                    {'sliderRadio',layerName, [{@obj.LaeryIO},{layerName}], {@obj.Transparancy, layerName}, '' });
                
            end
            
            obj.Plot;                                                       % Plot the new layer
            
        end
        
        
        % Hide a ploting layer that has been dismised --------------------------------------------------------------------------------------
        function obj = HideLayer(obj, layerName)
            
            % Make sure this function was not called out of progression
            allLayers = keys(obj.Layers);
            
            if ~any(ismember(keys(obj.Panels),'Layers'))
                disp('Hidding Layer befor Layers toolbar exists')
                return
            elseif ~any(ismember(allLayers,layerName))
                disp('Trying to Hide a none existant layer')
                return
            else
                
                % Remove the Layer from th plotting area
                obj.LayerOnOff(layerName) = false;
                obj.Plot;
                
                % Remove button group from Layers panel
                obj.GUIbuttons( 'Layers', 'remove', layerName);
                
            end
        end
        
        
        % Restor a hidden plotting layer
        function obj = RestoreLayer(obj, layerName)
            
            % Make sure this function was not called out of progression
            allLayers = keys(obj.Layers);
            
            if ~any(ismember(keys(obj.Panels),'Layers'))
                disp('Restoring Layer befor Layers toolbar exists')
                return
            elseif ~any(ismember(allLayers,layerName))
                disp('Trying to Restore a none existant layer')
                return
            else
                
                % Show the Layer in the plotting area
                obj.LayerOnOff(layerName) = true;
                obj.Plot;
                
                % Put the slider back button group from Layers panel
                obj.GUIbuttons( 'Layers','add', ...
                    {'sliderRadio',layerName, [{@obj.LaeryIO},{layerName}], {@obj.Transparancy, layerName}, '' });
            end
        end
        
        
        % Remove a plotting layer and its transparancy buttons ----------------------------------------------------------------------------
        function obj = DeleteLayer(obj, layer)
            
            allLayers = keys(obj.Layers);
            
            if any(ismember(allLayers,layer))               % make Sure the Layer exists
                remove(obj.Layers,           layer);        % Remove Plotting layer
                remove(obj.LayerTransparancy,layer);        % Remove layer's transparancy
                remove(obj.LayerOnOff,       layer);        % Remove layer IO
                
                obj.GUIbuttons( 'Layers', 'remove', layer); % Delete corresponding buttons and text
                
                obj.Plot;                                   % Replot data to clear the shape from the GUI display area
            end
            
        end
        
        
        % Make Figure -----------------------------------------------------------------------------------------------------------------------------
        function MakeFigure(obj, varargin)
            disp('Save Figure')
            
            fig = figure('units','normalized',...
                'Color',[1, 1, 1],...
                'NumberTitle','off',...
                'ToolBar','none' );
            
            
            copyobj(obj.Axes, fig);
            
            set(gca, 'Units', 'Normalized'); % First change to normalized units.
            set(gca, 'OuterPosition', [0, 0, 1, 1]);
            
        end
        
        
        % ____ Interactive functions ________________________________________________________________________________________________________________
        % Handel Interactive objects on the gui ---------------------------------------------------------------------------------------------------
        function pos = InteractiveLayer(obj, varargin)
            % varargin: array of varying size
            % * Could be empty or it could have latatude and longituded data to form the interactive boundry
            
            old_string = obj.Header.String;
            obj.Header.String = 'Doulbe Click Inside The Polygon to Accept and Close';                             % Prompt to select data
            
            % Disable all buttons so they dont interfear with the intteractive layer
            obj.DisalbeAllButtons;
            
            % Create interactive object and wait for it to be closed
            fcn = makeConstrainToRectFcn('impoly', [min(obj.Xgrid(:)), max(obj.Xgrid(:))], [min(obj.Ygrid(:)), max(obj.Ygrid(:))]);
            
            if nargin == 3
                iapg = impoly(obj.Axes, [varargin{1}, varargin{2}]);
            else
                iapg = impoly(obj.Axes);
            end
            
            setPositionConstraintFcn(iapg,fcn);                             % Create Interactive box with boarder constraints
            pos = wait(iapg);
            
            delete(iapg)                                                    % Delete the interactive Polygon
            
            obj.Header.String = old_string;                                 % But Back the origonal text
            
            % Enable the buttons again for continued user input
            obj.EnableAllButtons;
 
        end
               
        
        % Interactive Crop Box
        function [im, rect] = CropBox(obj, varargin)
            
            old_string = obj.Header.String;
            obj.Header.String = 'Doulbe Click Inside The Polygon to Accept and Close';                             % Prompt to select data
            
            % Disable all buttons so they dont interfear with the intteractive layer
            obj.DisalbeAllButtons;
            
            [im, rect] = imcrop(obj.Axes);
            
            obj.Header.String = old_string;                                 % But Back the origonal text
            
            % Enable the buttons again for continued user input
            obj.EnableAllButtons;
            
        end
        
        
        % Slidable Rectangle
        function pos = Slid_Rect(obj, varargin)
            % varargin: array of varying size
            % * Could be empty or it could have latatude and longituded data to form the interactive boundry
            
            old_string = obj.Header.String;
            obj.Header.String = 'Doulbe Click Inside The Polygon to Accept and Close';                             % Prompt to select data
            
            % Disable all buttons so they dont interfear with the intteractive layer
            obj.DisalbeAllButtons;
            
            % Create interactive object and wait for it to be closed
            fcn = makeConstrainToRectFcn('imrect', [min(obj.Xgrid(:)), max(obj.Xgrid(:))], [min(obj.Ygrid(:)), max(obj.Ygrid(:))]);
            
            if nargin == 3
                
                x = obj.Axes.XLim(2) / 2 - varargin{1}/2;
                y = obj.Axes.YLim(2) / 2 - varargin{2}/2;
                
                iapg = imrect(obj.Axes, [x, y varargin{1}, varargin{2}] );

                setResizable(iapg, false);                                  % Make the rectangle a fixed size
            
            else
                iapg = imrect(obj.Axes);
            end
            
            setPositionConstraintFcn(iapg, fcn);                             % Create Interactive box with boarder constraints
 
            pos = wait(iapg);
            
            delete(iapg)                                                    % Delete the interactive Polygon
            
            obj.Header.String = old_string;                                 % But Back the origonal text
            
            % Enable the buttons again for continued user input
            obj.EnableAllButtons;
            
        end
        
        % Mouse scrole wheel
        function obj = AddMouseWheel(obj, callBack, panel, min, max)
           obj.Figure.WindowScrollWheelFcn = {@obj.CountMouseWheel, callBack, panel, min, max};           
        end
        
        % Adjust Max Layer Transparancey
        function obj = AdjustAplahMax(obj, alpha_max)
            obj.AlphaMax = alpha_max;
        end
        
    end
    
    
    methods (Access = private)  %% Callback functions ===============================================================================================================
        
        % Layer Transparancy ----------------------------------------------------------------------------------------------------------------
        function obj = Transparancy(obj, source, event)
            obj.LayerTransparancy(source.Tag) = event.Source.Value * obj.AlphaMax;  % Get new layer transparancy and addjust from the normalized value
            obj.Plot;
        end
        
        
        % Turn layers onn and off ----------------------------------------------------------------------------------------------------------------
        function obj = LaeryIO(obj, source, event, ~)
            % Inputs: {'UIControl', 'Value change data', 'Layer Name' }
            
            obj.LayerOnOff(source.Tag) = event.Source.Value; % Get layer IO status
            obj.Plot;
            
        end
        
        
        % Finction to make sure you really want to close the GUI window
        function obj = CloseRequest(obj, varargin)                  
            selection = questdlg('Close This Figure?',...                   % Dialog box to ask user if they want to close the window
                'Close Request Function',...
                'Yes','No','Yes');
            
            switch selection
                case 'Yes'
                    close(obj.Figure)
                    
                case 'No'
                    return
            end
            
        end
        
        
        % Disable all buttons so they dont interfear with the intteractive layer
        function obj = DisalbeAllButtons(obj)
                        
            allPanels = keys(obj.Panels);                                   % Get Buton Panels In GUI Window
                      
            for i = 1:length(allPanels)                                     % Cycle through panels
  
                Panel = obj.Panels(allPanels{i});
                
                for j = 1:length(Panel.button)                              % Cycle throught button groups
                    Panel.button(j).DisableButtons;                         % Disabe Buttons
                end
            end
            
        end
        
        
        % Enable All buttons in gui window
        function obj = EnableAllButtons(obj)
            
            allPanels = keys(obj.Panels);                                   % Get Buton Panels In GUI Window

            for i = 1:length(allPanels)                                     % Cycle through panels
                
                Panel = obj.Panels(allPanels{i});
                
                for j = 1:length(Panel.button)                              % Cycle throught button groups
                    Panel.button(j).EnableButtons;                          % Disabe Buttons
                end
            end
        end
        
        
        % Count turns of the mouse scrole wheel
        function obj = CountMouseWheel(obj, ~, scroll, callBack, panel, min, max)
            
            obj.mouseCounter = obj.mouseCounter + scroll.VerticalScrollCount;   % Accumulate the scroll count
            
            if obj.mouseCounter < min                                       % Make sure the count is not less than the min
                obj.mouseCounter = min;
                return                                                      % Setting hasent changed, no need to exicute call back
            end
            
            if obj.mouseCounter > max                                       % Make sure the count does not exceed the max
                obj.mouseCounter = max;
                return                                                      % Setting hasent changed, no need to exicute call back
            end
            
            feval(callBack,panel);

                       
        end
        
    end
end

