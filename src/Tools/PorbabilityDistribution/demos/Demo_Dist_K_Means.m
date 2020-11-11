% Jakco Anderson
% Multivariate Gaussians K-means Clustering
% Sept 17, 2020
close all
clearvars
clc


n  = 100;                                                                   % Number of distribuitions
sz = 2;                                                                     % Distribuition Size ( dimentions )
k  = 10;                                                                    % Number of clusters



%% Create multivariant distributions

dist(n).mu  = rand(sz,1);                                                   % Pre alocate memory
dist(n).cov = zeros(sz);

for n = 1: n
    
    dist(n).mu  = rand(sz,1) * 30 - 15;
    dist(n).cov = zeros(sz);
    
    for a = 1:2:sz
        b = a+1;
        
        cov = rand(2,2) .* [2,.5;0.5,2] + eye(2) * rand;
        
        if rand < 0.6
            cov(2,1) = cov(1,2);
        else
            cov(2,1) = - cov(1,2);
            cov(1,2) = - cov(1,2);
        end
        
        dist(n).cov(a:b,a:b) = cov * cov / 2;
        
    end
end



%% --- Do clustering ----

indx = Dist_K_Means(dist, k, true);


