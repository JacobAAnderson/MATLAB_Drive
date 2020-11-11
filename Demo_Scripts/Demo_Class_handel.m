% Demo Class
% Jacob Anderson
% Sept 9, 2020


classdef Demo_Class_handel  < handle
      
    % Public Members accessable outside the class
    properties
        name
        gender
        mem1
        
    end
    
    
    % Private members only accessable inside the class
    properties (Access = private)
        mem2
        mem3
    end
    
    
    % Events to notify and listen to
    events
        NameChange
        MemberChange
        IdentityChange
    end
    
    
    % Public methods that can be accesses outside the class
    methods
        
        % Constructor --> Called when object is created
        function obj = Demo_Class_handel(name_, gender_)
            
            obj.name   = name_;
            obj.gender = gender_;
            
            obj.mem1 = 0;
            obj.mem2 = "Normal";                                % Access private properties from inside the class
            obj.mem3 = [];
            
            addlistener(obj, 'IdentityChange', @obj.Mixup);     % Add listender to an event
        end
        
        
        function obj = Change_Identity(obj, new_name, new_gender)
            obj.name = new_name;
            obj.gender = new_gender;
            
            if rand <= 0.5 
                notify(obj, 'IdentityChange')                               % Notify that an event has happend
            end
        end
        
        
        function obj = Change_Member(obj, new_member)
            obj.mem2 = new_member;
        end
        
        
        function obj = Change_Name(obj, new_name)
            obj.name = new_name;
            notify(obj,'NameChange')                                        % Notify that an event has happend
            
            obj = obj.Count_Changes;                                        % Call private method from inside the class
            
        end
        
        
        function obj = Change_Number(obj, number)
            
            arguments
                obj;
                number(:,:) {mustBeNumeric}
            end
            
            obj.mem3 = number; 
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
            notify(obj,'MemberChange')                                      % Notify that an event has happend
        end
        
        
        function obj = Mixup(obj, ~, ~)
            
            obj.name   = "I don't Know";
            obj.gender = "What would You Like it to be?";
            
            obj.mem1 = 96;
            obj.mem2 = "????";
            obj.mem3 = rand;
        end
        
        % Destructor --> called when object is deleted
        function delete(obj)
            
            disp("Pleas Saper Me!!")
            disp(obj)
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
        
        
        % This type of function only works with a handle object
        function obj = rando
            
            name_ = char( randi(100, 1,randi(20) ) );
            gend_ = char( randi(100, 1,randi(10) ) );
            
            obj = Demo_Class_handel( name_, gend_ );
            
            obj.mem1 = rand;
            obj.mem2 = "follow the rainbow!!";
            obj.mem3 = rand;
        end
    end
    
end