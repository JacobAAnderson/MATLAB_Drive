
close all
clear all
clc

speed = 2;
theta = 135;
t = 1;

sigma_Speed = 0.25;
sigma_theta = 1;

sigma_theta = sigma_theta * pi / 180;

sigma = zeros(2);

x = 0;
y = 0;

figure('name', 'Covariance')
plot(x,y, '+')
hold on

for ii = 1: 1

x = x + speed * t * cosd(theta);
y = y + speed * t * sind(theta);

sigS_2 = (sigma_Speed/speed) * (sigma_Speed/speed);

sig_x2 = sigS_2 + (sind(theta) * sigma_theta)^2;
sig_y2 = sigS_2 + (cosd(theta) * sigma_theta)^2;

sig_xy = -0.8*sind(2*theta) * sqrt(sig_x2) * sqrt(sig_y2);

theta = theta + 5;

mu = [x, y];

sigma = sigma + [sig_x2, sig_xy; sig_xy, sig_y2];

plot(x,y, '+')
plotErrorEllipse(mu, sigma, 90)

end

hold off