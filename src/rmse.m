



function err = rmse(X)

n = size(X,2);

X = X .* X;

err = sum(X,2)./n;

err = sqrt(err);

end