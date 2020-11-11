


function mse = MSError(A, B)

mse = [];

if any( size(A) ~= size(B))
    warning("Varialbe Must be the Same Size")
    return
end
    
dif = A - B;

dif = dif .* dif;

mse = nanmean(dif(:));


end