

close all
clear all
clc


for m =  0: 10 : 180

theta = m * pi/180;

x = [linspace(0,1,20); zeros(1,20)] + rand(2,20) * 0.1;


R = [ cos(theta), -sin(theta);
      sin(theta), cos(theta)];

X = R*x;


C = cov(X');


[V,D] = eig(C);

if V(2,1) > 0
    V(:,1) = -V(:,1);
    V(2,1) = -V(2,1);
end


x_ = mean( X(1,:) );
y_ = mean( X(2,:) );

u = [V(:,1); 0];
v = [V(:,2); 0];
c=cross(u,v);
nc=norm(c);
d=dot(u,v);
ThetaInDegrees = atan2d(nc,d)

theta = angleFromVectors(V(:,1), [0;-1])


l = [([0, V(1,2)] + x_ );
     ([0, V(2,2)] + y_ )];

s = [([0, V(1,1)] + x_ );
     ([0, V(2,1)] + y_ )];

figure('name', sprintf("M = %d", m))
plot(X(1,:), X(2,:),'.c')

hold on
plotErrorEllipse([x_, y_], C, 0.99)

% plot([0,V(1,1)], [0, V(2,1)], 'r' )
% plot([0,V(1,2) * sqrt(D(2,2)) ], [0, V(2,2) * sqrt(D(2,2)) ], 'm' )

plot(l(1,:), l(2,:), 'k' )
plot(s(1,:), s(2,:), 'm' )


hold off
axis equal
drawnow
end




function theta = angleFromVectors(a,b)


n1 = sqrt(a(1)*a(1) + a(2)*a(2));
n2 = sqrt(b(1)*b(1) + b(2)*b(2));

ab = a(1) * b(1) + b(2) * a(2);

a = ab/ (n1*n2);

theta = acos(a);

if abs( theta ) <= pi, return, end

if theta < 0
    warning( 'negative theta')
end

theta = 2*pi - theta;


end
