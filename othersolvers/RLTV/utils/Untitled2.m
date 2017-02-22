close all; clc
clear;
[X,Y] = meshgrid(-1:0.5:1);
Z = X.*Y;

distributionStep(X, Y, Z, 1, 4, 1, 'cube');