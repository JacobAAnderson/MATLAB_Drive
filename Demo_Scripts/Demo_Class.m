% Demo Class
% Jacob Anderson
% Sept 9, 2020


classdef Demo_Class
    
    % Public Members accessable outside the class
    properties
        name            % Class members that can take on values
        gender
        mem1
        
    end
    
    
    % Private members only accessable inside the class
    properties (Access = private)
        mem2
        mem3
    end

    
    % Public methods that can be accesses outside the class
    methods
        
        % Constructor
        function obj = Demo_Class(name_, gender_)
            obj.name   = name_;
            obj.gender = gender_;
            
            obj.mem1 = 0;
            obj.mem2 = "Normal";    % Access private properties from inside the class
            obj.mem3 = 0;
        end
        
        
        function obj = Change_Name(obj, new_name)
            obj.name = new_name;
            obj = obj.Count_Changes;    % Call private method from inside the class
        end
        
        
        function obj = Change_Number(obj, number)
            obj.mem3 = number; 
        end
        
        
        function obj = Change_Member(obj, new_member)
            obj.mem2 = new_member; 
        end
        
        
        function Show_Member(obj)
            disp(obj.mem2)
        end
        
        
        function Show_Number(obj)
            disp(obj.mem3)
        end
        
        
        % Only Save the data that is important
        function s = saveobj(obj)
            s.name   = obj.name;
            s.gender = obj.gender;
        end
        
        
    end
    
    
    % Methods that cannot be accessed from outside the class
    methods (Access = private)
        
        function obj = Count_Changes(obj)
            obj.mem1 = obj.mem1 + 1; 
        end
        
        
        function obj = Mixup(obj, ~, ~)
            
            obj.name   = "I don't Know";
            obj.gender = "What would You Like it to be?";
            
            obj.mem1 = 96;
            obj.mem2 = "????";
            obj.mem3 = rand;
        end
        
        
    end
    
    
    % Static methods that dont take any inputs
    methods (Static)
        
        % Load data that has been save
        function obj = loadobj(s)
            
            if isstruct(s)
                
                newObj =  Demo_Class(s.name, s.gender);
                
                newObj.mem1 = 0;
                newObj.mem2 = "Born Again";
                newObj.mem3 = 31;
                
                obj = newObj;
                
            elseif isa(s, 'Demo_Class')
                obj = s;
                
            else
                obj.name   = "I don't Know";
                obj.gender = "What would You Like it to be?";
                
                obj.mem1 = 96;
                obj.mem2 = "????";
                obj.mem3 = rand;
                
            end
        end
        
    end
    
end

