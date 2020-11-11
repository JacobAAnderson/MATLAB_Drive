% Image Window
% Jacob Anderson
% 9/3/2019

% Create a window inside an image without exceeding the boundaries


classdef SlidingWindow
    
    properties
        window_size
        sheet_size
    end
    
    methods
        
        function obj = SlidingWindow(window_size, sheet_size), warning on backtrace
            
            if nargin > 0
                obj = obj.SetDims(window_size, sheet_size);
            else
                obj.window_size = [];
                obj.sheet_size  = [];
            end
        end
        
        function obj = SetDims(obj,window_size, sheet_size)
            
            % Clear Previous Settings
            obj.window_size = [];
            obj.sheet_size  = [];
            
            if any( window_size > sheet_size)
                warning('Window Size Exceeds Sheet Size')
                return
            end
            
            % Assing new Settings
            obj.window_size = window_size;
            obj.sheet_size = sheet_size;
            
        end
        
        
        function rows = D1(obj, center)
            
            rows= [];
            
            if numel(obj.window_size) > 1 || numel(center) > 1, return, end
            
            h = round(obj.window_size/2);
            
            rows = center - h : center + h-1;
            
            if any(rows > obj.sheet_size)
                rows = rows - (rows(end) - obj.sheet_size + 1);
            end
            
            if any(rows < 1)
                rows = rows - rows(1) + 1;
            end
            
        end
        
        
        function [rows, cols] = D2(obj,center)
            
            m = obj.sheet_size(1);
            n = obj.sheet_size(2);
            
            if obj.window_size(1) > m || obj.window_size(2) > n
                rows= [];
                cols = [];
                return
            end
            
            r = center(1);
            c = center(2);
            
            h = round(obj.window_size/2);
            
            rows = r - h : r + h-1;
            cols = c - h : c + h-1;
            
            if any(rows > m)
                rows = rows - (rows(end) - m + 1);
            end
            
            if any(rows < 1)
                rows = rows - rows(1) + 1;
            end
            
            if any(cols > n)
                cols = cols - (cols(end) - n + 1);
            end
            
            if any(cols < 1)
                cols = cols - cols(1) + 1;
            end
            
        end
    end
end