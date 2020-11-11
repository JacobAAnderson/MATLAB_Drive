% Jacob Anderson
% Bhattacharyya Distance
% Sept 16, 2020



function distance = Bhattacharyya( dist1, dist2)

E = (dist1.cov + dist2.cov) /2;

m = dist1.mu - dist2.mu;

A = (m' * E^-1 * m);

B = log( det(E) / sqrt(det(dist1.cov) * det(dist2.cov)) );

distance = A/8 + B/2;

end
