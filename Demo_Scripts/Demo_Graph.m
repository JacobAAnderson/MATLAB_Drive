% Graph Demo
% Jacob Anderson
% 4/3/202

% Build a tree graph

close all
clear all
clc

cost = 0; cov = 0; action = false;

G = graph();
G = addnode(G,table(cost, cov, action));

leafs.indx = 1;
leafs.value = '0';
n = 1;

for t = 1:4
    for ii = size(leafs,2): -1:1    % Iterate through leafnoads
        
        p = leafs(ii).indx;
        
        for a = 1:3                     % Iterate through actions
        
            cov = rand;
            if cov < 0.7
                n = n+1;
                s = size(leafs,2) + 1;
                leafs(s).indx = n;
                leafs(s).value = string(n);
                
                if a == 2, cost = 1; action = 1;
                else, cost = 0; action = 0;
                end
                
                G = addnode(G,table(cost, cov, action));
                G = addedge(G,p,n);
               
            end
        end
        
        leafs([leafs.indx] == p) = [];
        
    end
    plot(G)
    drawnow
    
end


costs = G.Nodes.cost([leafs.indx]);
covs  = G.Nodes.cov([leafs.indx]);

in = costs == min(costs);
in = covs  == min(covs(in));

t = [leafs(in).indx];

P = shortestpath(G,1,t);
disp("path")
disp(P)

policy = G.Nodes.action(P);
disp("Policy")
disp(policy)



%% Constrianed tree with base nodes

close all
clear all
clc

% G = graph();
% G = addnode(G,16);
% G = addedge(G, [1,2,4,7,10], [2,4,7,10,13], 1);
% G = addedge(G, [3,5,6,8,9,11,12], [4,7,7,10,10,13,13], 1);
% G = addedge(G, [1,2,3,4,5,7,8,10,11,12], [3,5,6,8,9,11,12,14,15,16], 0);
% P = shortestpath(G,1,16);
% disp(P)
% plot(G)


leafs.indx = 1;                         % Node ID in Graph
leafs.value = '0';                      % Leaf Node value
leafs.count = 1;                        % Number of active leafnodes for prealocated memory

leafs(30).indx = 1;                     % Prealocat memory
leafs(30).indx = [];                    % Clear prealocat value


cov = 0; action = false;                % Data for the graph

G = graph();                            % Create Graph object
G = addnode(G,table(cov, action));      % Add root node

n = 1;                                  % Number of nodes in the graph

for t = 1:6
    
    s = [leafs.indx];
    
    action = false;
    
    for ii = leafs(1).count: -1:1                           % Iterate through leafnoads
        
        p = leafs(ii).indx;                                 % Parent Node being expanded
        
        cov = rand;
        if cov < 0.8
            
            n = n+1;
            G = addnode(G,table(cov, action));
            G = addedge(G,p,n,0);
            
            leafs(ii+1).indx = n;
            leafs(ii+1).value = string(n+1);
            
        else
            leafs(ii+1) = [];
            leafs(1).count = leafs(1).count - 1;            % Subtract a leaf node
        end
    end
    
    leafs(1).count = leafs(1).count + 1;                    % Add a leafnode
    
    cov = 0;
    action = true;
    
    G = addnode(G,table(cov, action));
    n = n+1;
    disp(n)
    
    l = repmat(n,size(s));
    
    G = addedge(G,s,l,1);
    
    leafs(ii).indx = n;
    leafs(ii).value = string(n);
    
    plot(G)
    drawnow
    
end

s = [leafs.indx];

p = zeros(size(s,2),t+1);

for ii = size(s,2): -1: 1    
    [p(ii,:),d(ii)] = shortestpath(G,1,s(ii));
end

covs  = G.Nodes.cov([leafs.indx]);
in = covs  == min(covs(d==min(d)));

P = p(in,:);
disp("path")
disp(P)

policy = G.Nodes.action(P);
disp("Policy")
disp(policy)





