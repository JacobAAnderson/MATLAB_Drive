% K-means distribuition clustering with Battacharia distance
% Jacob Anderson
% RDML, OSU, Corvallis, OR
% Oct 21, 2020




function indx = Dist_K_Means(dists, k, plotIO)

n = numel(dists);

% ----- 1) Randomly group the distributions into groups ------
C = uniqueIndex(n,k);                                                       % Pick k dirtibuitons to be the first centroids
centroids = dists(C);

centroid_last = centroids;                                                  % Keep track of last centroid to see how much it moved

zs = GetZer0s(centroids, k);                                                % Create zero centroid to calculate new centroids


C    = 1:k;                                                                 % Clusters
indx = NaN(n,1);                                                            % Prealocate cluster index

for t = 1:500
    
    % ----- 2) Group Distributions by centroid -----
    d = NaN(n,k);
    for a = 1:n
        for b = 1:k
%             d(a,b) = vecnorm( centroids(b).mu - dists(a).mu );                % < ----- !
            d(a,b) = Bhattacharyya( centroids(b), dists(a));                  % Calculate the distances between centroids and the distribuitions
        end
    end
    
    
    [~,indx] = min(d,[],2);                                                 % Get clustering index
   
%     PlotClusters(centroids, dists, indx, n, k)                              % Show new clusters
    
    % -- Deal with empty clusters ---                                       % Take the furthest pont from the largest cluster and assign it to a new cluster
    A = accumarray(indx,1,[k,1]);                                           % Number of nodes assigned to each cluster
    for a = find(A == 0)'                                                   % Check for empty clusters
            
            [~,i] = max(A);                                                 % Node with the most members
            
            mem = indx == i;                                                % Find members of that cluster
            
            cand = d(mem,i);                                                % Find the node that is furthest from the centroid
            cand = d(:,i) == max(cand);
            
            indx(cand) = a;                                                 % Put that node in the empty cluster
            
            A = accumarray(indx,1,[k,1]);                                   % Re-evaluate cluster size
            
    end
    
    
    % ----- 3) Calculate New Centroid -------
    
    centroids = zs;                                                         % Zero out the Centridos
    count = zeros(1,k);
    % Sum means and inverse covariance
    for a = 1:n
        
        in = indx(a) == C;
        centroids(in).mu  = centroids(in).mu  + dists(a).mu;                  % < ----- !
        centroids(in).cov = centroids(in).cov + dists(a).cov;
        count(in) = count(in) + 1;
        
%         centroid(in).mu  = centroid(in).mu  + sig * dist(a).mu;
%         sig = dist(a).cov^-1;
%         centroid(in).cov = centroid(in).cov + sig;
        
    end
    
    
    % Inver new covariances, normalize the new means and se if they moved
    shift = NaN(1,k);
    for a = 1:k
        
        centroids(a).mu  = centroids(a).mu  ./ count(a);                                % < ----- !
        centroids(a).cov = centroids(a).cov ./ count(a);
        
        shift(a) = Bhattacharyya( centroids(a), centroid_last(a));
        
%         shift(a) =  vecnorm( centroids(a).mu - centroid_last(a).mu);             % < ----- !
           
%         colors = GetColors;
%         plot(centroids(a).mu(1), centroids(a).mu(2), '+', 'color', colors(a,:) )
%         centroid(a).cov = centroid(a).cov^-1;
%         centroid(a).mu = centroid(a).cov * centroid(a).mu;
        
%         d(a) =  Bhattacharyya( centroid(a), centroid_last(a));
        
    end
    

    if all(shift < 0.0003)
%         fprintf('\n\nCentroids are stable after %d iterations\n\n', t)
        break
    end
    
    centroid_last = centroids;
    
end


% --- Make sure an empty cluster didn't slip by ---
% A = accumarray(indx,1,[k,1]);
% if any(A == 0), wrning('Empty Clusters Detected'), end


% --- Plot the Clusters ---
if nargin ==3 && plotIO
    
     PlotClusters(centroids, dists, indx, n, k)
                              
end

end



function C = uniqueIndex(n,k)

C = [];
while numel(C) < k
    C = randi(n,k,1);
    C = unique(C);
end

end


function zs = GetZer0s( centroids, k )

zs = centroids;

for t = 1:k
    zs(t).mu  = centroids(t).mu * 0;
    zs(t).cov = centroids(t).cov * 0;
end

end


function PlotClusters(centroids, dists, indx, n, k)

colors = GetColors;

figure('name', 'Distribution Clustering', 'numbertitle','off')

hold on

for ii = 1:n
     PlotDist( dists(ii), colors(indx(ii),:) )
%    plot(dists(ii).mu(1), dists(ii).mu(2),'o', 'color', colors(indx(ii),:) )
end

for ii = 1:k
    
    plot(centroids(ii).mu(1), centroids(ii).mu(2), '+', 'color', colors(ii,:) )
%     PlotDist( centroids(ii), [1, 0, 1] )
%     Plot_Ellipsoid(centroids(ii))
end

hold off
drawnow

end



function PlotDist(dist, color_)

p = 0.99;

mu  = dist.mu(1:2);
sig = dist.cov(1:2,1:2);


s = -2 * log(1 - p);

[V, D] = eig(sig * s);

t = linspace(0, 2 * pi);
a = (V * sqrt(D)) * [cos(t(:))'; sin(t(:))'];

a = real(a);

plot(a(1, :) + mu(1), a(2, :) + mu(2), 'color', color_)

end


function Plot_Ellipsoid(dist)

p = 0.99;

mu  = dist.mu(1:2);
sig = dist.cov(1:2,1:2);


s = -2 * log(1 - p);

[V, D] = eig(sig * s);


a = (V * sqrt(D)) * [1;1];


ellipsoid(mu(1), mu(2), 0, a(1), a(2), 1)

end



function c = GetColors(n)

if nargin == 1
    c = rand(n,3);
    
else
    c = [0                         0.447058823529412         0.741176470588235
        0.850980392156863         0.325490196078431         0.0980392156862745
        0.929411764705882         0.694117647058824         0.125490196078431
        0.494117647058824         0.184313725490196         0.556862745098039
        0.466666666666667         0.674509803921569         0.188235294117647
        0.635294117647059         0.0784313725490196        0.184313725490196
        1                         1                         0.0666666666666667
        1                         0.0745098039215686        0.650980392156863
        0                         0                         0
        0.0588235294117647        1                         1 ];
    
end

end