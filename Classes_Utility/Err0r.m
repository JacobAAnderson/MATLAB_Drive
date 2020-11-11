% Err0r Class to track errors
% Jacob Anderson
% RDML, OSU, Corvallis, OR.
% Nov 7, 2020


classdef Err0r
    
    properties
        base
        maxE
        minE
        threshold;
        idx
        
    end
    
    methods
        
        function obj = Err0r(sz, init)
            
            % Alocate memory
            obj.base(sz)      = 0;
            obj.maxE(sz)      = 0;
            obj.minE(sz)      = 0;
            obj.threshold(sz) = 0;
            
            
            % Assign Initial Value
            obj.base(1)      = init;
            obj.maxE(1)      = init;
            obj.minE(1)      = init;
            obj.threshold(1) = init;
            
            obj.idx = 1;
        end
        
        
        function obj = Add_Datum(obj, base_, max_, min_, threshold_)
            
            obj.idx = obj.idx + 1;
        
            obj.base(obj.idx)      = base_;
            obj.maxE(obj.idx)      = max_;
            obj.minE(obj.idx)      = min_;
            obj.threshold(obj.idx) = threshold_;
            
        end
        
        
        function Plot(obj)
            
            t = obj.idx;
            
            p0 = plot(obj.base(1:t));
            hold on
            p1 = plot(obj.maxE(1:t));
            p2 = plot(obj.minE(1:t));
            p3 = plot(obj.threshold(1:t));
            hold off
            legend([p0,p1,p2,p3], {'Base','max','Min','thresh'},'Location','northwest')
            title("Expected Error")
            xlabel("Time steps [s]")
            ylabel("Error [m]")
            
        end
        
        
        
    end
    
end




