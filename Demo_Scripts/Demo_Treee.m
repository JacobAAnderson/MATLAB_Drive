% Demo Treee --> Tree structure graph
% Jacob Anderson
% 2/18/2020

clear all
close all
clc



node.Action = 'S';
node.XY = [0,0];
node.Err = 0;

tree = Treee('Tree Planning');
tree = tree.RootNode(node);

actions = ['S', 'L', 'R', 'F', 'B'];                     % Posible Actions
move    = [   0,0;   -1,0;     1,0;       0,1;   0,10];                     % Coresponding Motions

T = 10;

% Grow Tree
for t = 1:T
    
    leafs = tree.LeafNodes;         % Indecies of the leaf nodes
    
    for ii = 1:numel(leafs)
        
        for jj = 1:randi(3)
            
            xy  = tree.Values(leafs(ii), 'XY');
            err = tree.Values(leafs(ii), 'Err');
            
            act = randi(numel(actions));
            
            node.Action = actions(act);
            node.XY = xy + move(act,:);
            node.Err = err + 0.1;
            
            [tree, ~] = tree.AddNode(leafs(ii), node);
            
        end
        
    end
end


tree.Plot('Err')


%%
clc 
clear all
close all

load('tree.mat')

tree.Plot('Error')




%%
clc

leafs = tree.LeafNodes;

times = tree.LeafValues('Time');
error = tree.LeafValues('Error');

idx = times == max(times);

leafs = leafs(idx);
error = error(idx);

idx = error == min(error);

idx = leafs(idx);

policy = tree.Search(idx(1), 'Action');


disp(policy)

