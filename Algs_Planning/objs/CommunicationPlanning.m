% Communication Planning
% Jacob Anderson
% 2/11/2020


classdef CommunicationPlanning
    
    properties
        name
        stats
        vehicles
        
    end
    
    methods
        
        function obj = CommunicationPlanning(name)
            obj.name = name;
            obj.vehicles.count = 0;
            %             obj.Tree = struct('State', [], 'Cost', [], 'Parent', []);
        end
        
        
        % Add vehicles to the planning algorithm
        function obj = Add_Vehicle(obj, pf, path, name_, r, maxspeed, speedNoise, compassNoise)
            
            % Check if the info is packed into a cell from being used inside the function call
            if nargin == 2
                
                path         = pf{2}; 
                name_        = pf{3}; 
                r            = pf{4}; 
                maxspeed     = pf{5};
                speedNoise   = pf{6}; 
                compassNoise = pf{7};
                pf           = pf{1};
                
            end
            
            
            
            ii = obj.vehicles(1).count + 1;
            
            obj.vehicles(ii).PF         = pf;
            obj.vehicles(ii).path       = path;
            obj.vehicles(ii).name       = name_;
            obj.vehicles(ii).maxSpeed   = maxspeed;
            obj.vehicles(ii).wp_radious = r;
            obj.vehicles(ii).wp_idx     = 1;
            obj.vehicles(ii).CompassErr = compassNoise;
            obj.vehicles(ii).SpeedErr   = speedNoise;
            obj.vehicles( 1).count      = ii;
            
            obj.stats(ii).action.cov = [];
            
        end
        
        
        % Plan Communications - Origonal
        function policy = Plan0(obj, T, dt, cov_max, plotIO)
            
            if nargin < 5, plotIO = false; end
            
%             disp('Communications Planning')
            tic
            
            leafNodes(2*T).B = obj.vehicles;
            leafNodes(1).B = obj.vehicles;
            leafNodes(1).C = 0;
            leafNodes(1).E = 0;
            leafNodes(1).P = 0;
            leafNodes(1).I = 1;
            
            node.Action = false;
            node.Parent = 0;
            node.Time   = 0;
            node.Cost   = 0;
            node.Error  = 0;
            
            
            if plotIO
                figure('name', 'Simulaiton Tree', 'numbertitle','off')
                plot(0,0,'.g')
                hold on
            end
            
            
            tree = Treee('Comms Planning');                                 % Create Tree
            tree = tree.RootNode(node);                                     % Add Root Node to the tree
            
            
            % Find Parent Node for the Base Line --> Parent with lowst cost and uncertainty
            C = 0;                                                          % Node with lowest Cost
            P = 1;
            last = 0;
            
            count = 1;
            err_b = 100;
            
            num_leafs(T) = NaN;
            for t = 1:T
                
                num_leafs(t) = count;
                
                % Without Communications
                for ii = count : -1: 1
                    [B, err] = Expand(leafNodes(ii).B, dt);
                    
                    if err < cov_max * err_b
                        
                        node.Action = false;
                        node.Cost   = leafNodes(ii).C;
                        node.Time   = t;
                        node.Error  = err;
                        node.Parent = leafNodes(ii).I;
                        
                        [tree, idx] = tree.AddNode(node);
                        
                        leafNodes(ii+1).B = B;
                        leafNodes(ii+1).C = leafNodes(ii).C;
                        leafNodes(ii+1).E = err;
                        leafNodes(ii+1).P = leafNodes(ii).I;
                        leafNodes(ii+1).I = idx;
                        
                        if plotIO, plot(t, err, '.r'), end
                        
                    else
                        count = count - 1;
                    end
                    
                end
                count = count + 1;
                
                % With Communications
                [leafNodes(1).B, err_b, done] = BaseLine(leafNodes(1).B, dt);
                
                node.Action = true;
                node.Cost   = C+1;
                node.Error  = err_b;
                node.Parent = P;
                node.Time   = t;
                
                if plotIO, plot(t, err_b, '.g'), end
                [tree, idx] = tree.AddNode(node);
                
                leafNodes(1).C = C+1;
                leafNodes(1).E = err_b;
                leafNodes(1).P = P;
                leafNodes(1).I = idx;
                
                if done, break, end                                         % End Simulation if the vehicle has acchived its goal
                
                % Find Parent Node for the Base Line --> Parent with lowst cost and uncertainty
                I  = [leafNodes.I] < last;
                leafNodes(I) = [];
                
                C_ = [leafNodes.C];                                         % Node with lowest Cost
                P_ = C_ == min(C_);
                min_err = [leafNodes(P_).E];                                % Node with lowest cost and uncertainty
                P = [leafNodes.E] == min(min_err);
                P = leafNodes(P).I;
                C = min(C_);
                
                last = max([leafNodes.I]);
                
            end
            
%             save(sprintf('%s_%s_Comms_Policy_Tree.mat', TODAY, obj.Name),'tree')
            
            leafs = tree.LeafNodes;
            
            times = tree.LeafValues('Time');
            error = tree.LeafValues('Error');
            
            idx = times == max(times);
            
            leafs = leafs(idx);
            error = error(idx);
            
            idx = error == min(error);
            idx = leafs(idx);
            
            policy = tree.Search(idx(1), 'Action');
            
            fprintf('\nCommunications Planning Finnished\nTime: %f [min]\n',toc/60)
            fprintf('Number of Communications: %d\nExpected Error: %f\n',sum(policy), min(error))
            fprintf('Min Ending Error: %f\n\n', leafNodes(1).E)
            
            save('Num_leafs0.mat', 'num_leafs')
        end
        
        
        % Plan Communications - sdt clustering
        function policy = Plan1(obj, T, dt, cov_max, plotIO)
            
            k  = 10;        % Number of Clusters 
            k_ = 50;        % Clustering threshold
            
            if nargin < 5, plotIO = false; end
            
            % --- Simulation Queue ----
            leafs.indx   = 1;                                               % Node ID in Graph
            leafs.parent = 0;                                               % Parent Node
            leafs.value  = obj.vehicles;                                    % Leaf Node value
            leafs.count  = 1;                                               % Number of active leafnodes for prealocated memory
            leafs.cov    = 6;                                               % State uncertainty
            
            leafs(2*T).indx = 1;                                            % Prealocat memory
            leafs(2*T).indx = [];                                           % Clear prealocat value
            % ------------------------
            
            n   = 1;                                                        % Current node index
            cov = 6;                                                        % Current state uncertainty
            action = false; time = 1;                                       % Current communications acction

            G = digraph();                                                    % Instanciate the communications graph
            G = addnode(G,table(time,cov, action));                              % Add root node for t = 1

            if plotIO
                figure('name', 'Simulaiton Tree', 'numbertitle','off')
                plot(G)
            end
            
            err_b = 100;                                                    % Base Node Error
            
            num_leafs(T) = NaN;                                             % Keep track of the number of leaf nodes
            
            error_base(T)   = 0;
            error_max(T)    = 0;
            error_min(T)    = 0;
            error_thresh(T) = 0;
            error_base(1)   = 3;
            error_max(1)    = 3;
            error_min(1)    = 3;
            error_thresh(1) = 3;
            
%             tic
            for t = 2:T
                
                tic
                
                time = t;
                
                num_leafs(t) = leafs(1).count;
                
                % --- Simulate Without Communications ---
                s = [leafs.indx];                                           % List of leaf nodes to become the parents of the next base node
                action = false;                                             % No comms --> Action = false
                
                for ii = leafs(1).count : -1: 1                             % Iterate through leaf nodes
                   [B, cov, ~] = Expand(leafs(ii).value, dt);                % Forward simulate node without communication
%                    [B, cov] = Expand_all(leafs(ii).value, dt);             % Forward simulate node without communication and all otehr robots communicating
                    
                    % Compare the resulting localization uncertainty
                    if cov < cov_max * err_b || cov <= 25 % -> The New Simulation state is Accepted                            
                        
                        n = n+1;                                            % New node index
                        p = leafs(ii).indx;                                 % Parent Node index
                        G = addnode(G,table(time, cov, action));            % Add node to graph
                        G = addedge(G,p,n,0);                               % Add edge between parent node and new node
                        
                        leafs(ii+1).indx = n;                               % Add resulting simulation to the queue
                        leafs(ii+1).parent = leafs(ii).indx;
                        leafs(ii+1).value = B;
                        leafs(ii+1).cov = cov;
                        
                        if plotIO, plot(G), end
                        
                    else                                                    % -> The New Simulation State Is Not Accepted
                        leafs(ii+1) = [];                                   % Subtract a leaf node
                        leafs(1).count = leafs(1).count - 1;
                    end
                    
                end
                
                % --- Simulate With Communications ---
                n = n+1;                                                    % New node index
                [leafs(1).value, err_b, ~, done] = BaseLine(leafs(1).value, dt);   % Forward simulate base node with communication
                leafs(1).indx = n;                                          % Add to base of leafnodes
                leafs(1).count = leafs(1).count + 1;
                leafs(1).cov = err_b;
                
                l = repmat(n,size(s));                                      % Set up edges for the graph
                action = true;                                              % Node information
                cov = err_b;
                
                G = addnode(G,table(time, cov, action));                    % Add new node to the graph
                G = addedge(G,s,l,1);                                       % Add edges from all of the previous leafnodes to this base node
                
                error_base(t)   = err_b;
                error_max(t)    = max([leafs.cov]);
                error_min(t)    = min([leafs.cov]);
                error_thresh(t) = err_b * cov_max;
                
                % --- Merge Dublicate States ---
                if leafs(1).count >= k_
                   
                    idx = [0; kmeans([leafs(2:end).cov]',k)];               % Bin the leaf nodes by uncertainty and Skip the base node
                    
                    in = zeros(k+1,1);                                      % Index of the nodes to be saved

                    for bb = 1:k                                            % Iterate through the bins to rewire the graph
                        
                        id = idx == bb;                                     % Find all nodes in the current bin
                        id(1) = false;                                      % Keep base node safe
                        
                        p = [leafs(id).parent];                             % Parents of the binned nodes
                        c = [leafs(id).indx];                               % Graph index value of the binned nodes 
                        
                        G = rmedge(G,p,c);                                  % Remove the edges between parens and curent leaf nodes
                         
                        c = c(randi(numel(c)));                             % Pick leaf node to save
                        in(bb) = c;
                        
                        c = repmat(c,size(p));                              % Create edges from parent nodes to the saved node
                        G = addedge(G,p,c,0);
                        
                        
                    end
                    
                    in(bb+1) = leafs(1).indx;                               % Add base node into the final group of nodes to keep                   
                    
                    out = setdiff([leafs.indx],in);                         % Find nodes to illiminate
                    [~,aa] = ismember(out,[leafs.indx]);
                    leafs(aa) = [];                                         % Remove nodes
                    
                    leafs(1).count = k+1;                                   % Reset node count
                end
                
                
                if plotIO, plot(G), drawnow, end
                
                
                Write2File(sprintf('%s_PlanningAlg_stats.txt', obj.name), '%d, %d, %f\n', [t, leafs(1).count, toc] )
                
                fprintf('%s, T: %d, N: %d, Loop Time: %f\n',obj.name, t, leafs(1).count, toc)
                
                % --- Check For Goal ---
                if done, break, end                                         % End Simulation if the vehicle has acchived its goal
            end
            
            % --- Find Communications Policy ----
%             covs = G.Nodes.cov([leafs.indx]);
%             in = covs  == min(covs);
%             t = [leafs(in).indx];
%             P = shortestpath(G,1,t);                                        % Get path from optimal leafnode

            
            for ii = leafs(1).count: -1: 1 
                try
                [P(ii,:), d(ii)] = shortestpath(G,1,leafs(ii).indx);            % Get path from optimal leafnode
                catch
                    disp('')
                    pause
                end
            end
            
            covs = G.Nodes.cov([leafs.indx]);
            covs_min = covs(d == min(d));
            
            P = P(min(covs_min) == covs,:); 
            
            if size(P,2) ~= t
                disp('???')
                pause
            end
            
            policy = G.Nodes.action(P(1,:));                                     % Get Communications Policy
            
            figure('name',obj.name)
%             subplot(3,1,1)
%             plot(G.Nodes.time(P(1,:)))
            subplot(2,1,1)
            plot(policy)
            hold on
            
            p_ = double(policy);
            p_(p_==0) = NaN;
            plot(p_, '*')
            title( sprintf('%s Communication Policy',obj.name))
            xlabel('Time Steps [30 sec]')
            ylabel('Communicate [true/false]')
            
            
            subplot(2,1,2)
            p0 = plot(error_base(1:t));
            hold on
            p1 = plot(error_max(1:t));
            p2 = plot(error_min(1:t));
            p3 = plot(error_thresh(1:t));
            hold off
            title('Expected Error from Simulations')
            xlabel('Time Steps [30 sec]')
            ylabel('Error [m]')
            
            legend([p0,p1,p2,p3], {'Base','max','Min','thresh'},'Location','northwest')
            drawnow
            
            fprintf('\nCommunications Planning Finnished for %s\n', obj.name) %Time: %f [min]\n',toc/60)
            fprintf('Number of Communications: %d\n', sum(policy))
%            fprintf('Expected Error: %f\n\n', covs(in))
            
%             save('Num_leafs.mat', 'num_leafs')
        end
        
        
        % Plan Communications - distribuition clustering
        function policy = Plan2(obj, T, dt, cov_max, plotIO)
            
            n    = 1;               % Current node index
            k    = 10;              % Number of Clusters
            k_   = 50;              % Clustering threshold
            cov  = 6;               % Current state uncertainty 
            time = 1;               % Start Time
            err = 100;              % Base Node Error
            action = false;         % Current communications acction
            
            if nargin < 5, plotIO = false; end
            
            leafs = Queue(2*T, {'index','parent','sim','cov','dist'});     % Simulation Queue
            leafs = leafs.Add_Node(1, 0, obj.vehicles, 6, []);              % Add the first node
            
            error = Err0r(T,3);                                             % Keep track of results
            
            
            G = digraph();                                                  % Instanciate the communications graph
            G = addnode(G,table(time,cov, action));                         % Add root node for t = 1
            
            
            if plotIO   % Debugging outputs
                figure('name', 'Simulaiton Tree', 'numbertitle','off')
                plot(G)
            end
            
            
            for t = 2:T
                
                tic
                
                time = t;
                
                s = leafs.List('index');                                    % List of leaf nodes to become the parents of the next base node
                
                
                % --- Forward Simulate Without Communications -------------
                action = false;                                             % No comms --> Action = false
                
                for ii = leafs.count : -1: 1                                % Iterate through leaf nodes
                    
                    [B, cov, dist] = Expand(leafs.Value(ii, 'sim'), dt);    % Forward simulate node without communication
                    
                    % --> The New Simulation state is Accepted
                    if cov < cov_max * err || cov <= 25                     
                        
                        n = n+1;                                            % New node index
                        p_idx = leafs.Value(ii, 'index');                       % Parent Node index
                        
                        G = addnode(G,table(time, cov, action));            % Add node to graph
                        G = addedge(G,p_idx,n,0);                               % Add edge between parent node and new node
                        
                        leafs = leafs.Replace_Node(ii+1, n, p_idx, B, cov, dist);   % Add new node to the queue
                        
                        if plotIO,  plot(G), end                             % Debugging outputs, plot the graph
                    
                        
                    % --> The New Simulation State Is Not Accepted
                    else, leafs = leafs.Remove_Node(ii+1);                  % Subtract a leaf node
                    end
                    
                end
                
                % --- Forward Simulate With Communications ----------------
                [B, err, dist, done] = BaseLine(leafs.Value(1, 'sim'), dt); % Forward simulate base node with communication
                
                n = n+1;  
                
                leafs = leafs.Replace_Node(1, n, [], B, err, dist);         % Update the base node
                                                        
                leafs.count = leafs.count + 1;
                
                
                l = repmat(n,size(s));                                      % Set up edges for the graph
                action = true;                                              % Node information
                cov = err;
                
                G = addnode(G,table(time, cov, action));                    % Add new node to the graph
                G = addedge(G,s,l,1);                                       % Add edges from all of the previous leafnodes to this base node
                
                
                % Keep track of performance
                error = error.Add_Datum(err, max(leafs.List('cov')), min(leafs.List('cov')), err * cov_max);
                

                % --- Merge Dublicate States ------------------------------
                if leafs.count >= k_
                    
                    % 1) Cluster nodes by particle filter distribuition
                    idx = [0; leafs.Cluster_Nodes(k, true)];                % Perform Node clustering and add the base node to the clusters
                    
                    in = zeros(k+1,1);                                      % Zero out the index of nodes to keep (k+1 because we save the base)
                    
                    
                    % 2) Iterate throught the clusters, select a node from each cluster to save and discard the rest
                    for bb = 1:k
                        
                        % 2.a) Remove the edges between parens and curent leaf nodes
                        id = idx == bb;                                     % Find all nodes in the current bin
                        id(1) = false;                                      % Keep base node safe
                        
                        p_idx = leafs.List('parent',id);                    % Parents of the binned nodes
                        c_idx = leafs.List('index', id);                    % Graph index value of the binned nodes
                        
                        G = rmedge(G,p_idx,c_idx);                          % Remove the edges between parens and curent leaf nodes
                        
                        
                        % 2.b) Pick node to save and create edges from parent nodes to the saved node
                        i = randi( numel(c_idx) );                          % Pick leaf node to save
                        c_ = c_idx(i);                             
                        in(bb) = c_;                                        % Add it to the list of saved nodes
                        
                        c_ = repmat(c_,size(p_idx));                        % Create edges from parent nodes to the saved node
                        G = addedge(G,p_idx,c_,0);
                        
                        
                        % 2.c) Remove the discarded nodes from the graph --> Problem with node ids changing in the graph
%                         c_idx(i) = [];                                      % Take the saved node out the cluster index
%                         G = rmnode(G,c_idx);                                % Remove the discarded nodes from the graph
                      
                            
                    end
                    
                    % 3) Remove the discarded nodes from the queue
                    in(end) = leafs.Value(1, 'index');                      % Add base node into the final group of nodes to keep
                     
                    out = setdiff( leafs.List('index'), in );               % Find nodes to illiminate
                    
                    [~,aa] = ismember(out, leafs.List('index') );
                    
                    leafs = leafs.Remove_Node(aa);                          % Remove nodes
                    
                    leafs.count = k+1;                                      % Reset node count
                end
                
                
                if plotIO, plot(G), drawnow, end
                
                leafs.growth(t) = leafs(1).count;                           % Keep track of howmany simulations are present
                
%                 Write2File(sprintf('%s_PlanningAlg_stats.txt', obj.name), '%d, %d, %f\n', [t, leafs(1).count, toc] )
                
%                 fprintf('%s, T: %d, N: %d, Loop Time: %f\n',obj.name, t, leafs(1).count, toc)
                
                
                % --- Check For Goal ---
                if done, break, end                                         % End Simulation if the vehicle has acchived its goal
                
            end
            
            
            % --- Find Communications Policy ------------------------------
            for ii = leafs(1).count: -1: 1
                try
                    [P(ii,:), d(ii)] = shortestpath(G,1,leafs.Value(ii, 'index') );    % Get path from optimal leafnode
                catch E
                    warning(E.message)
                    pause
                end
            end
            
            covs = G.Nodes.cov(leafs.List('index'));
            covs_min = covs(d == min(d));
            
            P = P(min(covs_min) == covs,:);
            
            if size(P,2) ~= t
                disp('???')
                pause
            end
            
            policy = G.Nodes.action(P(1,:));                                     % Get Communications Policy
            
            figure('name',obj.name)
            subplot(3,1,1)
            plot(G.Nodes.time(P(1,:)))
            title("time Steps")
            subplot(3,1,2)
            plot(policy)
            title("Communication Ploicy")
            subplot(3,1,3)
            error.Plot;
            drawnow
            
            fprintf('\nCommunications Planning Finnished for %s\n', obj.name) %Time: %f [min]\n',toc/60)
            fprintf('Number of Communications: %d\n', sum(policy))
%             fprintf('Expected Error: %f\n\n', covs(in))
%             
%             save('Num_leafs.mat', 'num_leafs')
        end
        
        
        % Plan Communications - Monte Carlo Tree Search
        function policy = Plan_MCTS(obj, T, dt, cov_max, plotIO)
            
            if nargin < 5, plotIO = false; end
            
            % --- Simulation Queue ----
            leafs.indx   = 1;                                               % Node ID in Graph
            leafs.parent = 0;                                               % Parent Node
            leafs.value  = obj.vehicles;                                    % Leaf Node value
            leafs.count  = 1;                                               % Number of active leafnodes for prealocated memory
            leafs.cov    = 6;                                               % State uncertainty
            
            leafs(2*T).indx = 1;                                            % Prealocat memory
            leafs(2*T).indx = [];                                           % Clear prealocat value
            % ------------------------
            
            n   = 1;                                                        % Current node index
            cov = 6;                                                        % Current state uncertainty
            action = false; time = 1;                                       % Current communications acction
            
            G = digraph();                                                  % Instanciate the communications graph
            G = addnode(G,table(time,cov, action));                         % Add root node for t = 1
            
            if plotIO
                figure('name', 'Simulaiton Tree', 'numbertitle','off')
                plot(G)
            end
            
            num_leafs(T) = NaN;                                             % Keep track of the number of leaf nodes
            
            error_base(T)   = 0;
            error_max(T)    = 0;
            error_min(T)    = 0;
            error_thresh(T) = 0;
            error_base(1)   = 3;
            error_max(1)    = 3;
            error_min(1)    = 3;
            error_thresh(1) = 3;
            
            %             tic
            for t = 2:T
                
                time = t;
                num_leafs(t) = leafs(1).count;
                
                % --- Simulate Without Communications ---
                s = [leafs.indx];                                           % List of leaf nodes to become the parents of the next base node
                action = false;                                             % No comms --> Action = false
                
                for ii = leafs(1).count : -1: 1                             % Iterate through leaf nodes
                    [B, cov] = Expand(leafs(ii).value, dt);                % Forward simulate node without communication
                    %                    [B, cov] = Expand_all(leafs(ii).value, dt);             % Forward simulate node without communication and all otehr robots communicating
                    
                    % Compare the resulting localization uncertainty
                    if (cov < cov_max * err_b && cov >= err_b) || cov <= 20 % -> The New Simulation state is Accepted
                        
                        n = n+1;                                            % New node index
                        p = leafs(ii).indx;                                 % Parent Node index
                        G = addnode(G,table(time, cov, action));            % Add node to graph
                        G = addedge(G,p,n,0);                               % Add edge between parent node and new node
                        
                        leafs(ii+1).indx = n;                               % Add resulting simulation to the queue
                        leafs(ii+1).parent = leafs(ii).indx;
                        leafs(ii+1).value = B;
                        leafs(ii+1).cov = cov;
                        
                        if plotIO, plot(G), end
                        
                    else                                                    % -> The New Simulation State Is Not Accepted
                        leafs(ii+1) = [];                                   % Subtract a leaf node
                        leafs(1).count = leafs(1).count - 1;
                    end
                    
                end
                
                % --- Simulate With Communications ---
                n = n+1;                                                    % New node index
                [leafs(1).value, err_b, done] = BaseLine(leafs(1).value, dt);   % Forward simulate base node with communication
                leafs(1).indx = n;                                          % Add to base of leafnodes
                leafs(1).count = leafs(1).count + 1;
                leafs(1).cov = err_b;
                
                l = repmat(n,size(s));                                      % Set up edges for the graph
                action = true;                                              % Node information
                cov = err_b;
                
                G = addnode(G,table(time, cov, action));                    % Add new node to the graph
                G = addedge(G,s,l,1);                                       % Add edges from all of the previous leafnodes to this base node
                
                error_base(t)   = err_b;
                error_max(t)    = max([leafs.cov]);
                error_min(t)    = min([leafs.cov]);
                error_thresh(t) = err_b * cov_max;
                
                if plotIO, plot(G), drawnow, end
                
                % --- Check For Goal ---
                if done, break, end                                         % End Simulation if the vehicle has acchived its goal
            end
            
            % --- Find Communications Policy ----
            %             covs = G.Nodes.cov([leafs.indx]);
            %             in = covs  == min(covs);
            %             t = [leafs(in).indx];
            %             P = shortestpath(G,1,t);                                        % Get path from optimal leafnode
            
            
            for ii = leafs(1).count: -1: 1
                try
                    [P(ii,:), d(ii)] = shortestpath(G,1,leafs(ii).indx);            % Get path from optimal leafnode
                catch
                    disp('')
                    pause
                end
            end
            
            covs = G.Nodes.cov([leafs.indx]);
            covs_min = covs(d == min(d));
            
            P = P(min(covs_min) == covs,:);
            
            if size(P,2) ~= t
                disp('???')
                pause
            end
            
            policy = G.Nodes.action(P(1,:));                                     % Get Communications Policy
            
            figure('name','time')
            subplot(3,1,1)
            plot(G.Nodes.time(P(1,:)))
            subplot(3,1,2)
            plot(policy)
            subplot(3,1,3)
            p0 = plot(error_base(1:t));
            hold on
            p1 = plot(error_max(1:t));
            p2 = plot(error_min(1:t));
            p3 = plot(error_thresh(1:t));
            hold off
            legend([p0,p1,p2,p3], 'Base','max','Min','thresh')
            drawnow
            
            fprintf('\nCommunications Planning Finnished\n') %Time: %f [min]\n',toc/60)
            fprintf('Number of Communications: %d\n', sum(policy))
            %            fprintf('Expected Error: %f\n\n', covs(in))
            
            %             save('Num_leafs.mat', 'num_leafs')
        end
        
        
        % On the fly planning
        function [obj, comms] = Update_Belief(obj, dt, self)
            
            for ii = 1: obj.vehicles(1).count
                
                if strcmp(obj.vehicles(ii).name, 'self'), continue, end
                
                pf0 = obj.vehicles(ii).PF;
                
                [speed, heading, obj.vehicles(ii).wp_idx, ~] = WaypointFollow(pf0.X, obj.vehicles(ii).path, obj.vehicles(ii).wp_radious, obj.vehicles(ii).wp_idx);
                
                speed = speed * obj.vehicles(ii).maxSpeed + random(obj.vehicles(ii).SpeedErr,1);            % Conver Speed from normalized value to vehicle speed
                
                heading = CompassAngle(heading, 'rad', 'deg') + random(obj.vehicles(ii).CompassErr,1);      % Conver heading to a compass angleHeading
                
                dr = DeadReckoning(pf0.X(1:2)', speed, heading, dt, 'deg', 'compass');
                
                measurment = pf0.Map.Depth(dr) + normrnd( 0, 0.06);
                
                
                % ___ Update with out Comms ______________________________
                pf1 = pf0.Update(speed, dt, heading, measurment);
                
                % ___ Update with comms __________________________________
                dist = vecnorm( pf0.X(1:2) - self.X(1:2));
                
                acoms = {'acoms', dist + random('Normal', 0, 0.5), self.X(1:2)', self.cov(1:2,1:2)};
                
                pf2 = pf0.Update(speed, dt, heading, measurment, acoms);
                
                
                if trace(pf1.cov(1:2,1:2)) > 1.25 * trace(pf2.cov(1:2,1:2))
                    comms = true;
                    obj.vehicles(ii).PF = pf2;
                else
                    comms = false;
                    obj.vehicles(ii).PF = pf1;
                end
                
            end
            
            
        end
        
        
        
        % Make Sure that our simulated position is within the bounds of our actual position
        function [obj, comms] = Check_Belief(obj, dt, self)
            
            
            % ___ Forward Simulate Self --> Constriants ___________________
            names = {obj.vehicles.name};
            ii = find(strcmp(names, 'self'));
            pf0 = obj.vehicles(ii).PF;
            
            [speed, heading, obj.vehicles(ii).wp_idx, ~] = WaypointFollow(pf0.X, obj.vehicles(ii).path, obj.vehicles(ii).wp_radious, obj.vehicles(ii).wp_idx);
            
            speed = speed * obj.vehicles(ii).maxSpeed  + random(obj.vehicles(ii).SpeedErr,1);           % Conver Speed from normalized value to vehicle speed
            
            heading = CompassAngle(heading, 'rad', 'deg') + random(obj.vehicles(ii).CompassErr,1);      % Conver heading to a compass angleHeading
            
            dr = DeadReckoning(pf0.X(1:2)', speed, heading, dt, 'deg', 'compass');
            
            measurment = pf0.Map.Depth(dr) + normrnd( 0, 0.06);
            
            pf0 = pf0.Update(speed, dt, heading, measurment);
            
            
            if vecnorm(self.X(1:2)' - pf0.X(1:2)') > trace(self.cov(1:2,1:2))
                comms = true;
                obj.vehicles(ii).PF = self;
            else
                comms = false;
                obj.vehicles(ii).PF = pf0;
            end
            
            
            
            
        end
        
        
    end
end





% Base Line Step, Vehicle communicates
function [vehicles, err, dist, done] = BaseLine(vehicles, dt)

names = {vehicles.name};
me = find(strcmp(names, 'self'));

if me ~= 1                                                                  % Make Sure that I am first in line
    self_ = vehicles(me);
    vehicles(me) = vehicles(1);
    vehicles(1) = self_;
end

done = false(vehicles(1).count,1);
err = 0;

dist.mu  = zeros( vehicles(1).count *2, 1);
dist.cov = zeros( vehicles(1).count *2 );

c = vehicles(1).count;  

for ii = [1:c; 1:2:c*2]
    
    if ii(1) == 1
        [vehicles(1), done(1)] = Simulate(vehicles(1), dt);
        X   = vehicles(1).PF.X;
        Cov = vehicles(1).PF.cov;
%         err = trace(vehicles(ii).PF.cov);
    else
        [vehicles( ii(1) ), done( ii(1) )] = Simulate(vehicles( ii(1) ), dt, X, Cov);
%         if trace(err) < trace(vehicles(ii).PF.cov), err = vehicles(ii).PF.cov; end
    end
    
    set = ii(2): ii(2)+1;
    dist.mu(set)      = vehicles( ii(1) ).PF.X(1:2);
    dist.cov(set,set) = vehicles( ii(1) ).PF.cov(1:2,1:2);
    
    err = err + trace(vehicles( ii(1) ).PF.cov);
end

% err  = trace(err);
done = any(done);

end


% Expand Node with out communication
function [vehicles, err, dist] = Expand(vehicles, dt)

err = 0;

dist.mu  = zeros( vehicles(1).count *2, 1);
dist.cov = zeros( vehicles(1).count *2 );

c = vehicles(1).count;  

for ii = [1:c; 1:2:c*2]
    
    [vehicles( ii(1) ), ~] = Simulate(vehicles( ii(1) ), dt);
    
    set = ii(2): ii(2)+1;
    dist.mu(set)      = vehicles( ii(1) ).PF.X(1:2);
    dist.cov(set,set) = vehicles( ii(1) ).PF.cov(1:2,1:2);
    
    
    err = err + trace(vehicles( ii(1) ).PF.cov);

%     err_ = trace(vehicles(ii).PF.cov);
% 
%     if err_ > err, err = err_; end
    
%      if ~strcmp(vehicles(ii).name, 'self') && trace(err) < trace(vehicles(ii).PF.cov)
%         err = vehicles(ii).PF.cov;
%     end
end
end


% Simulate A Vehicle following waypoints
function [Vehicle, done] = Simulate(Vehicle, dt, X, Cov)

pf = Vehicle.PF;

[speed, heading, Vehicle.wp_idx, done] = WaypointFollow(pf.X, Vehicle.path, Vehicle.wp_radious, Vehicle.wp_idx);

speed = speed * Vehicle.maxSpeed  + random(Vehicle.SpeedErr,1);                % Conver Speed from normalized value to vehicle speed

heading = CompassAngle(heading, 'rad', 'deg') + random(Vehicle.CompassErr,1);  % Conver heading to a compass angleHeading

dr = DeadReckoning(pf.X(1:2)', speed, heading, dt, 'deg', 'compass');

measurment = pf.Map.Depth(dr) + normrnd( 0, 0.06);                          % Get measurment and add sensor noise

if nargin < 3
    pf = pf.Update(speed, dt, heading, measurment);
else
    dist = vecnorm( pf.X(1:2) - X(1:2)) + random('Normal', 0, 0.5);
    
%     acoms = {'acoms', dist + random('Normal', 0, 0.5), X(1:2)', Cov(1:2,1:2)};
    
    pf = pf.Acoms(dist, X(1:2)', Cov(1:2,1:2));
    pf = pf.Update(speed, dt, heading, measurment);
    
end

Vehicle.PF = pf;

end




%% Old Stuff  


% function [vehicles, err] = Expand_all(vehicles, dt)
% 
% err_b = inf;
% 
% for a = 0: vehicles(1).count
%     
%     if a == 1, continue, end
%     
%     [vehicles_, err_, ~] = Evaluate(vehicles, dt, a);
%     
%     if err_ < err_b 
%         vehicles = vehicles_;
%         err = err_;
%     end 
%     
% end
% 
% end


% function [vehicles, err, done] = Evaluate(vehicles, dt, action)
% 
% done = false(vehicles(1).count,1);
% err = 0;
% 
% if action > 0
%    X = vehicles(action).PF.X;
%    Cov = vehicles(action).PF.cov;
% end
% 
% for ii = 1: vehicles(1).count
%     
%     if action == 0 || action == ii
%         [vehicles(ii), done(ii)] = Simulate(vehicles(ii), dt);
%     else
%         [vehicles(ii), done(ii)] = Simulate(vehicles(ii), dt, X, Cov);
%     end
%     
%     err_ = trace(vehicles(ii).PF.cov);
% 
%     if err_ > err, err = err_; end
%     
% end
% 
% done = any(done);
% 
% end




