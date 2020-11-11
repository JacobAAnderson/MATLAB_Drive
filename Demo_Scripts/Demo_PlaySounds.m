%% Play sounds
clear all
close all
clc

load chirp
sound(y, Fs)
plot(y)
%%
load gong
sound(y, Fs)
plot(y)

%%
load handel
sound(y, Fs)
plot(y)
%%
load laughter
sound(y, Fs)
plot(y)
%%
load mtlb
sound(mtlb, Fs)
plot(y)

%%
load splat
sound(y, Fs)
plot(y)

%%
load train
sound(y, Fs)
plot(y)


