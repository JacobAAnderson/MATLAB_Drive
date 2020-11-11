% Tree
% Jacob Anderson
% 2/18/2020

classdef Treee
    
    properties
        Name
        Nodes
        Index
    end
    
    methods
        
        
        function obj = Treee(name), warning on backtrace
            obj.Name = name;
            obj.Index = 1;
        end
        
        
        function obj = RootNode(obj, node)
            obj.Nodes = node;
            obj.Nodes.Parent = 0;
            obj.Nodes(10000) = node;                                        % Preallocate memory
        end
        
        
        function [obj, idx] = AddNode(obj, node)
            idx = obj.Index + 1;
            obj.Nodes(idx) = node;
            obj.Index = idx;
        end
        
        
        function leafNodes = LeafNodes(obj)
            parents = unique([obj.Nodes.Parent]);
            nodes = 1: obj.Index;
            leafNodes = setdiff( nodes, parents);
        end
        
        
        function values = LeafValues(obj, domain)
            idx = obj.LeafNodes;
            values = obj.Values(idx,domain);
        end
        
        
        function values = Values(obj, idx, domain)
            
            values(size(idx,1), size(obj.Nodes(1).(domain),2)) = 0;
            
            for ii = 1: numel(idx)
                values(ii,:) = obj.Nodes(idx(ii)).(domain);
            end
            
        end
        
        
        function policy = Search(obj, start, domain)
            
            policy = [];
            
            if ~isfield(obj.Nodes(1), domain), warning('This Field Is Not Present'), return, end
            if start > obj.Index, warning("Starting Node is outside of this tree's domain"),return, end
            
            idx = obj.BackSearch(start);
            policy = idx;
            
            for ii = 1: numel(idx)
                policy(ii) = obj.Nodes(idx(ii)).(domain);
            end
            
            policy = fliplr(policy);
            
        end
        
        
        
        function Plot(obj, param)
            
            figure
            
            parents = [obj.Nodes.Parent];
            
            y = obj.Nodes(1).(param);
            
            plot(0,y,'*g')
            hold on
            
            for ii = 1: 20
                
                idx = find(parents == ii);
                
                for jj = 1: numel(idx)
                    
                    x1 = obj.Nodes(idx(jj)).Time;
                    y1 = obj.Nodes(idx(jj)).(param);
                    
                    p = obj.Nodes(idx(jj)).Parent;
                    
                    x2 = obj.Nodes(p).Time;
                    y2 = obj.Nodes(p).(param);
                    
                    plot(x1,y1,'*r')
                    line([x1,x2],[y1,y2], 'color', 'c')
                    
                end
                
            end
        end
                
            
            
        
        
    end
    
    methods (Access = private)
        
        function path = BackSearch(obj, idx)
            
            path = zeros(1,1000);                                           % Preallocate memory
            ii = 1;
            
            while idx
                path(ii) = idx;
                idx = obj.Nodes(idx).Parent;
                ii = ii+1;
            end
            
            path(ii:end) = [];                                               % Get Rid of extra space
            
        end
        
    end
end
