
function plotErrorEllipse(mu, Sigma, p, color)

    s = -2 * log(1 - p);

    [V, D] = eig(Sigma * s);

    t = linspace(0, 2 * pi);
    a = (V * sqrt(D)) * [cos(t(:))'; sin(t(:))'];

    a = real(a);
    
    if nargin > 3
        plot(a(1, :) + mu(1), a(2, :) + mu(2), color)
    else
        plot(a(1, :) + mu(1), a(2, :) + mu(2));
    end
    
end