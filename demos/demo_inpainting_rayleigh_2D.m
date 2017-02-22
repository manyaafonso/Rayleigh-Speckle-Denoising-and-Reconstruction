close all;
clear;

addpath(genpath('../src'))
addpath(genpath('../data'))
addpath(genpath('../utils'))
addpath(genpath('../othersolvers'))


x0 = imresize( double(rgb2gray(imread('lena512color.tiff'))), 0.5 );
[M,N] = size(x0);

x = abs(x0)/max(x0(:));

%%%% generate mask
fractionmissing = 0.5;
q = rand(M,N);
mask = double(q>fractionmissing);

%%%% observation
y0 = mask.*x;
noise0 = sqrt(raylrnd(ones(M,N)));
y = noise0.*y0;

figure, imagesc(x0), colormap gray, axis image, axis off;

figure, imagesc(y0), colormap gray, axis image, axis off;
figure, imagesc(y), colormap gray, axis image, axis off;


mu1 = 50;
mu2 = 50;
lambda = 2;

opts = struct('maxiters',200,'chambolleit',5,'inner_iters',5,'stopcriterion',2,'tol', 1e-2,'verbose',0,'x_true',x);
[clean, obj, times, mses] = rayleighReconstruct2D_v2(y,mask,lambda,mu1,mu2,opts);
fprintf('Iters: %g, CPU time: %g seconds, MSE: %g\n', length(times)-1,times(end),mses(end))
x_hat = clean;
mae_admm = sum(abs(x(:)-x_hat(:)))/numel(x);
fprintf('ADMM\nIters: %g, CPU time: %g seconds, MAE: %g\n', length(times)-1,times(end),mae_admm)

figure, imagesc(x_hat), colormap gray, axis image, axis off;



alfa = lambda; % regularization prior
nIter = 50;

cx = repmat([1:M]',1,N);
cy = repmat([1:N],M,1);
cx = cx.*mask;
cy = cy.*mask;

cg = y;

t0 = cputime;
[mlr,a,b,xn,yn] = maxLh_interp2D(cx, cy, cg, max(M,N));
sig = usDenoising(mlr,alfa,nIter);
x_rltv = sig./max(sig(:));

tfinal = cputime-t0;

mae_rltv = sum(abs(x(:)-x_rltv(:)))/numel(x);
fprintf('RLTV\nIters: %g, CPU time: %g seconds, MAE: %g\n', nIter,tfinal,mae_rltv)

figure, imagesc(sig), colormap gray, axis image, axis off;

