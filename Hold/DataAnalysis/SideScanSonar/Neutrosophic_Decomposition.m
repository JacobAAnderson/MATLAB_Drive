% Neutrosophic Decomposition
% Jacob Anderson
% May 17, 2019

function [T, F , I] = Neutrosophic_Decomposition(image)

% Tuening Parameters ---------------------------------------------
alpha_min = 0.01;
alpha_max = 0.1;


% Primamry Neutrosophic Decomposition -------------------------------
mf = fspecial('average',5);

T = imfilter(image,mf,'replicate');
I = abs( image - T);

T = ( T - min(T(:)) ) / ( max(T(:)) - min(T(:)) );
I = ( I - min(I(:)) ) / ( max(I(:)) - min(I(:)) );
F = 1 - T;


% decrease subset indetermancy --------------------------------------
T_prim = 2*T.*T;
T_prim(T > 0.5) = 1-2*(1-T(T > 0.5)).*(1-T(T > 0.5));

F_prim = 2*F.*F;
F_prim(F > 0.5) = 1-2*(1-F(F > 0.5)).*(1-F(F > 0.5));


[pixelCounts, ~ ] = imhist(image);                                          
pdf = pixelCounts / numel(image);                                           % Computer probability density function

H_i = - nansum( pdf .* log(pdf) );
H_max = log( numel(image) );

alpha = alpha_min + (alpha_max - alpha_min) * H_i / H_max;
beta = 1-alpha;

T(I >= beta) = T_prim(I>= beta);
F(I >= beta) = F_prim(I>= beta);

end

