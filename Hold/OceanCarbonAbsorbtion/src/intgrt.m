function int=intgrt(z, x)

% M-file function to integrate  x wrt z using trapeziodal rule.

[r,c]=size(z);
int=[0];
for i=1:(r-1)
    int = int + (z(i+1)-z(i))*(x(i+1)+x(i))/2;
end

