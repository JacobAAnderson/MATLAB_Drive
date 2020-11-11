% Simulation Queue
% Jacob Anderson
% RDML, OSU, Corvallis, OR.
% Nov 6, 2020



classdef Queue
    
    properties
        nodes
        datum
        count
        growth
    end
    
    
    methods
        
        % Instanciate the Queue
        function obj = Queue(T, properties)
            
            for p = properties
                
                obj.nodes(T).(p{1}) = [];                                      % Create node structure and prealocate memory
                
            end
           
            obj.count  = 0;                     % Number of active leafnodes for prealocated memory
            obj.growth = NaN(T/2,1);                     % Number of nodes at each time step
            
        end
        
 
        % Add a node to the queue
        function obj = Add_Node(obj, varargin)
            
            obj.count  = obj.count + 1;                                     % Increase the node count
            
            idx = obj.count;
            
            for ii = [fields(obj.nodes(idx))' ; varargin]
                
                obj.nodes(idx).(ii{1}) = ii{2};
            end
            
        end
               
        
        % Replace a node to the queue
        function obj = Replace_Node(obj, idx, varargin)
            
            for ii = [fields(obj.nodes(idx))' ; varargin]
                
                obj.nodes(idx).(ii{1}) = ii{2};
            end
            
        end
        
        
        % Remove a node from the queue
        function obj = Remove_Node(obj, idx)
            
            obj.nodes(idx) = [];
            
            obj.count  = obj.count - numel(idx);
            
        end
        
        
        % Get Node value
        function v = Value(obj, idx, param), v = obj.nodes(idx).(param);  end
        
        
        % List of inxeces
        function s = List(obj, param, id) 
            
            if nargin == 3
                s = [obj.nodes(id).(param)]; 
            else
                s = [obj.nodes.(param)];
            end
            
        end
        
        
        % Node Clustering
        function idx = Cluster_Nodes(obj, k, dispIO)
            
            idx = Dist_K_Means([obj.nodes(2:end).dist], k, dispIO);       % Bin the leaf nodes by uncertainty
           
        end
        
        
        
    end
    
    
    
    
end