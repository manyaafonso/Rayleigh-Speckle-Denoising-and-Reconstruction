close all;
clear all;

addpath(genpath('../src'))
addpath(genpath('../data'))
addpath(genpath('../utils'))
addpath(genpath('../othersolvers'))

x0 = phantom3d(64);
x = x0/max(x0(:));
noise0 = sqrt(raylrnd(ones(size(x0))));

figure, view3Dmodified(x0)

y = x.*noise0;
y = y/max(y(:));

lambda = 1;
mu = 2.5;
opts = struct('maxiters',200,'chambolleit',5,'inner_iters',4,'stopcriterion',1,'tol', 1e-1,'verbose',0,'x_true',x);

[x_admm, obj_admm, times_admm, f1, f2] = rayleighDenoise3D_v1(y, lambda, mu, opts);
fprintf('Proposed\nIterations: %d, total CPU time: %g\n', length(times_admm)-1, times_admm(end))

figure, view3Dmodified(x_admm)
