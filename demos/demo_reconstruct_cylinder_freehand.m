close all;
clear all;

addpath(genpath('../src'))
addpath(genpath('../data'))
addpath(genpath('../utils'))
addpath(genpath('../othersolvers'))


n = 60;
noPlanes = 30;

sigma = 0.01;

% %%% comment these lines and uncomment load to save time
% v0 = createObject([n,n,n],0,1,'tube',0);
% v = raylrnd(v0);
% 
% [cx,cy,cz,cg] = sampleObject(v, n, noPlanes);
% close all;
% 
% save 'cylinder3dFH.mat' n v0 v noPlanes cx cy cz cg
% %%%% end of block to comment

load cylinder3dFH.mat;

%%%%%%% Rearrange samples and sample grid as a masking operation

[M,N,L] = size(v0);

xd = cx;
yd = cy;
zd = cz;
pd = cg;


xd = xd(:); yd = yd(:); zd = zd(:); pd = pd(:); pd(find(isnan(pd))) = 0;

xd = xd(:); yd = yd(:); zd = zd(:); 
pd = pd(:); pd(find(isnan(pd))) = 0;

ind = (xd<=1).*(xd>=-1).*(yd<=1).*(yd>=-1).*(zd<=1).*(zd>=-1);

xd = xd(ind==1);
yd = yd(ind==1);
zd = zd(ind==1);
pd = pd(ind==1);


% the coordinates of each plane are re-scaled to [-1,1] according to the 'rv' dimensions 
xn   = (2*xd - min(xd) - max(xd))./(max(xd)-min(xd));
yn   = (2*yd - min(yd) - max(yd))./(max(yd)-min(yd));
zn   = (2*zd - min(zd) - max(zd))./(max(zd)-min(zd));

% find closest nodes for each observed pixel after obs.
tx = ((xn - min(xn))*(n-1))./(max(xn)-min(xn))+1 ; 
tx = round(tx);

ty = ((yn - min(yn))*(n-1))./(max(yn)-min(yn))+1 ; 
ty = round(ty);

tz = ((zn - min(zn))*(n-1))./(max(zn)-min(zn))+1 ; 
tz = round(tz);


observed_voxels = [tx, ty, tz];


observed = zeros(n,n,n);
pixel_count = zeros(n,n,n);

tic
for i = 1:length(pd(:))
    
    if pixel_count(tx(i),ty(i),tz(i)) ==0
        
        observed(tx(i),ty(i),tz(i)) = pd(i);
        pixel_count(tx(i),ty(i),tz(i)) = 1;
    else
        observed(tx(i),ty(i),tz(i)) = (pd(i)+pixel_count(tx(i),ty(i),tz(i))*observed(tx(i),ty(i),tz(i)))/( pixel_count(tx(i),ty(i),tz(i))+1 );
        pixel_count(tx(i),ty(i),tz(i)) = pixel_count(tx(i),ty(i),tz(i)) + 1;
    end
    
end

mask = double(pixel_count>0);
y = observed;

%%%%%% reconstruct


mu1 = 50;
mu2 = 50;
alpha = 80;

opts = struct('maxiters',100,'chambolleit',5,'inner_iters',5,'stopcriterion',2,'tol', 1e-2,'verbose',0);
[clean, obj, times, mses] = rayleighReconstruct3D_v2(y, mask, alpha, mu1, mu2, opts);

mae = sum(abs(v0(:)-clean(:)))/numel(v0);

fprintf('Proposed: Iters: %g, CPU time: %g seconds, MAE: %g\n', length(times)-1,times(end),mae)

figure, plot(times,obj,'LineWidth',2.2), set(gca,'FontName','Times'), set(gca,'FontSize',14), xlabel('seconds');

clean = reshape(clean,M,N,L);

figure, view3Dmodified(v0)
figure, view3Dmodified(y)
figure, view3Dmodified(clean)


alfa = 1; % regularization prior
nIter = 100;
t0 = cputime;
[mlr,a,b,xn,yn,zn] = maxLh_interp(xd, yd, zd, pd, n);
[sig, iters_rltv] = usDenoising_v3(mlr,alfa,nIter,1e-4);
sig = sig./max(sig(:));

tfinal = cputime-t0;
mse_rltv = norm(sig(:)-v0(:),2)^2/numel(sig);
mae_rltv = sum(abs(sig(:)-v0(:)))/numel(sig);
fprintf('RLTV 3D: %d iters, %g seconds, MAE: %g\n', iters_rltv, tfinal, mae_rltv)
figure, view3Dmodified(sig)

