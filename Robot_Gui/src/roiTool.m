classdef roiTool < handle
    properties
        InteractiveObj
        pos
        tag
        
        funcs = {@(h) impoly(h),     'impoly';              % Interactive Polygon
                 @(h) imellipse(h),  'imellipse';           % Interactive ellipse
                 @(h) imfreehand(h), 'imfreehand';          % interactive free hand drawing
                 @(h) imline(h),     'imline';              % Interactive line
                 @(h) imdistline(h), 'imline';              % Interactive line tool with a distance measurement
                 @(h) impoint(h),    'impoint';             % Interactive points
                 @(h) imrect(h),     'imrect' };            % Interactive Rectangle
    end
    
    methods
        
        % Create Interactive objects on the gui ---------------------------------------------------------------------------------------------------
        function obj = roiTool(Axes, varargin)
            
            varargin = varargin{1};         % Extract input arguments from their cell
            
            obj.InteractiveObj = impoly(Axes, [varargin{1}, varargin{2}]);                                      % Creat interacitve polygon
            
            fcn = makeConstrainToRectFcn('impoly', [min(LON(:)), max(LON(:))], [min(LAT(:)), max(LAT(:))]);     % Constrain the polygon the image
            setPositionConstraintFcn(obj.InteractiveObj,fcn);
            
            addNewPositionCallback(obj.InteractiveObj,@obj.NewPos);     % Callback function to update the pose of the polygon
            
            if nargin == 5
                addNewPositionCallback(obj.InteractiveObj,varargin{3});
            end
            
            setColor(obj.InteractiveObj,color)    % Change the color of the polygone to make it distenct
            
            source = dbstack(2);        % Get calling function and make a note of it for futue refrance
            obj.tag = source(1).name;
            
        end
        
        
        % Call back function that updates the pose of the Interactive Object ---------------------------------------------------------------------
        function obj = NewPos(obj,varargin)
            obj.pos = getPosition(obj.InteractiveObj);
        end
        
        % Delete Interactive Object --------------------------------------------------------------------------------------------------------------
        function obj = DelteObj(obj,varargin)
            
            delete(obj.InteractiveObj)
            
        end
        
    end
end

