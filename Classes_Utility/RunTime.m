

classdef RunTime < handle
    
    properties % (Access = private)
        s
    end
    
    methods (Static)
        
        function obj = RunTime, warning on backtrace  
            tic
            obj.s = sprintf('%.1f',toc);
            fprintf("Elapsed Time: %s", obj.s)
        end
        
    end
    
    
    methods
        function Step(obj)
            fprintf(repmat('\b',1,numel(obj.s)));
            
            obj.s = sprintf('%.1f',toc);
            fprintf(obj.s)
            
        end
        
        
        function time = Time(~)
            time = toc;
        end
        
    end
end