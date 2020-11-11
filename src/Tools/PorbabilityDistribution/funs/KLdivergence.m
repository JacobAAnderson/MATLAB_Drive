% Jacob Anderson
% KL Divergence between two multivariate Gaussians
% Sept 17, 2020


function distance = KLdivergence(dist1, dist2)

A = log( det(dist2.cov)/ det(dist1.cov) );

B = trace( dist2.cov^-1 * dist1.cov );

C = (dist2.mu - dist1.mu)' * dist2.cov^-1 * (dist2.mu - dist1.mu);

distance = (A + B + C) / 2;

end