% Distribuiton Clustering
% Jacob Anderson
% Sept 16, 2020

close all
clear all
clc

format long g


dist1.mu  = rand(2,1);
dist1.cov = rand(2,2);

dist2.mu  = rand(2,1);
dist2.cov = rand(2,2);



dist3.mu  = rand(4,1);
dist3.cov = zeros(4);
dist3.cov(1:2,1:2) = rand(2,2);
dist3.cov(3:4,3:4) = rand(2,2);


         
dist4.mu  = rand(4,1);
dist4.cov = zeros(4);
dist4.cov(1:2,1:2) = rand(2,2);
dist4.cov(3:4,3:4) = rand(2,2);         


% Create large multivariant distributions
sz = 20;

dist5.mu  = rand(sz,1);
dist5.cov = zeros(sz);

dist6.mu  = rand(sz,1);
dist6.cov = zeros(sz);

for a = 1:2:sz
    b = a+1;
    dist5.cov(a:b,a:b) = rand(2,2);
    dist6.cov(a:b,a:b) = rand(2,2);
    
end



t = 10000;
         
times = zeros(t,6);
         
for i = 1:t
    
    tic
    d1 = Bhattacharyya( dist1, dist2);
    times(i,1) = toc;
    
    tic
    da = KLdivergence( dist1, dist2);
    times(i,2) = toc;
    
    % disp(' ')
    
    tic
    d2 = Bhattacharyya( dist3, dist4);
    times(i,3) = toc;
    
    tic
    db = KLdivergence( dist3, dist4);
    times(i,4) = toc;
    
    % disp(' ')
    
    tic
    d3 = Bhattacharyya( dist5, dist6);
    times(i,5) = toc;
    
    tic
    dc = KLdivergence( dist5, dist6);
    times(i,6) = toc;
    
end

disp( mean(times,1)')

