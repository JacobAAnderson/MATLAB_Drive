% Demo Kalman Filter
% Jacob Anderson
% 2/7/2020

clear all
close all
clc

a = [1, 0;
     0, 1];
b = [1, 0;
     0, 1];
 
c = [1, 0;
     0, 1];


kf = Kalman_Filter(a, b, c);
kf.Q = [0.5, 0.01; 0.01, 0.5];
kf.R = [ 0.3, 0; 0, 0.3];

clear a b c

X = [1;1];
sig = [0.1, 0; 0, 0.1];
u = [2;3];

z = [4;5];

[X, sig] = kf.Update(X,sig,u, z)
