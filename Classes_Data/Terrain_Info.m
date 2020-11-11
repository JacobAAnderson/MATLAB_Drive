

classdef Terrain_Info
    
    properties
        Map
        Window_Size
    end
    
    methods
        
        function obj = Terrain_Info(map), warning on backtrace 
            obj.Map = map;
            obj.Map.Name = 'Terrian Info';
        end
        
        
        function obj = Eval_MapInfo(obj, win_sz)
            
            % ___ Entropy, STD, Ave Elevation _____________________________
            obj.Window_Size = win_sz;
            
            sz = size(obj.Map.Map);
            
            sw = SlidingWindow(win_sz, sz); 
            
            entopy    = zeros(sz);
            stand_dev = zeros(sz);
            average   = zeros(sz);
            
            for row = 1 : sz(1)
                for col = 1 : sz(2)
                    
                    [rows, cols] = sw.D2([row, col]);
                    
                    if any( isnan( obj.Map.Map(rows, cols)))
                        entopy(   row, col) = NaN;
                        stand_dev(row, col) = NaN;
                        average(  row, col) = NaN;
                        
                    else
                        e = obj.Map.Map(rows, cols);
                        
                        p = e ./ nansum(e(:));
                        p = p .* log(p);
                        
                        entopy(   row, col) = -sum(p(:));
                        stand_dev(row, col) =  std(e(:));
                        average(  row, col) = mean(e(:));
                        
                    end
                end
            end
            
            entopy(   isnan(obj.Map.Map)) = NaN;
            stand_dev(isnan(obj.Map.Map)) = NaN;
            average(  isnan(obj.Map.Map)) = NaN;
            
            obj.Map = obj.Map.Add_Layer("Entropy", entopy);
            obj.Map = obj.Map.Add_Layer("STD", stand_dev);
            obj.Map = obj.Map.Add_Layer("Ave", average);
            
            % ___ Group by percent of points at theat elevation ___________
%             idx = kmeans(obj.Map.Map(:),10);
%             idx = reshape(idx, size(obj.Map.Map));
%             percent = idx;
%             tot = sum(~isnan(obj.Map.Map));     % Get the number of non nan values in the map
%             
%             for ii = 0:10
%                 s = sum(idx == ii);
%                 percent(idx == ii) = s/tot;
%             end
%             
%             obj.Map = obj.Map.Add_Layer("Cluster", percent);
            
            
            % ___ Second derivative of the terrain ________________________
            [FX,FY] = gradient(obj.Map.Map);
            
            FXX = gradient(FX);
            [~,FYY] = gradient(FY);
            
            grad =  (FXX + FYY);
            
            obj.Map = obj.Map.Add_Layer("Grad", grad);
            
            
        end
        
        
    end
    
end