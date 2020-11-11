% Datetime ismember function
% Jacob Anderson
% RDML, OSU, Corvallis OR.
% Oct 2, 2020

% Compare Date times while avoiding round off error


classdef DateTime_Funs
    
    properties
       errorThrehold 
    end
    
    methods
        
        function obj = DateTime_Funs(err_threshold)
            
            if nargin
                obj.errorThrehold = err_threshold;
            else
                obj.errorThrehold = 0.00000001;
            end
        end
        
        
        function tf = Isequal(obj, set1, set2)
            
            tf = ~ (abs(set1 - set2) > obj.errorThrehold);
            
        end
        
        
        function [Lia,Locb] = Ismember(obj, dateTime_A, dateTime_B)
            
            Lia  = false(size(dateTime_A));
            Locb = zeros(size(dateTime_A));
            
            for ii = 1: numel(dateTime_A)
                for jj = 1: numel(dateTime_B)
                    
                    if abs(dateTime_A(ii) - dateTime_B(jj)) < obj.errorThrehold
                        Lia(ii) = true;
                        Locb(ii) = jj;
                        continue
                    end
                end
            end
            
            Locb(Locb == 0) = [];
            
        end
        
        
        function list = Unique(obj, list)
            
            for ii = numel(list): -1: 1
                
                for jj = ii-1: -1: 1
                    
                    if obj.Isequal( list(ii), list(jj) )
                        list(ii) = [];
                        break
                    end
                end
            end
        end
        
        
        function tf = SetError(obj, set1, set2)
            
            if ~ all( obj.Isequal( set1, set2))
                dt = set1 - set2;
                dt.Format = 's';
                
                if nargout == 0
                    warning("\n\tDate-time Set error\n")
                    disp(dt)
                    
                else
                    tf = true;
                end
            
                
            elseif nargout > 0
                tf = false;
            end
            
        end
        
    end
    
    
    
end