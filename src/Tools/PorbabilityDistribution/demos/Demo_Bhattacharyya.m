% Demo Bhattacharyya Distance
% Jacob Anderson
% RDML, OSU, Corvallis, OR
% Oct 21, 2020




dist1.mu = 0.5;
dist1.cov = 1.6;

dist2.mu = 1.3;
dist2.cov = 0.7;

distance = Bhattacharyya( dist1, dist2);
disp(distance)


dist1.mu = [0.5; 1.2];
dist1.cov = [1.6, 0.6; 0.6, 2.3];

dist2.mu = [1.3; 2.2];
dist2.cov = [0.7, 0.9; 0.9, 0.8];

distance = Bhattacharyya( dist1, dist2);
disp(distance)




