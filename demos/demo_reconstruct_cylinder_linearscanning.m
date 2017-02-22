close all; 
clear all;

addpath(genpath('../src'))
addpath(genpath('../data'))
addpath(genpath('../utils'))
addpath(genpath('../othersolvers'))

% Create the 'original' ultrasound data volume
n = 60;
global v;
v0 = createObject([n,n,n],0,1,'tube',0);

v = raylrnd(v0);

% Simulate acquisition of 2D images
noPlanes = 30; % user-specified number of planes to acquire

mask = ones(n,n,n);
cx = zeros(n,n,noPlanes);
cy = zeros(n,n,noPlanes);
cz = zeros(n,n,noPlanes);

for i = 1:noPlanes
    
    mask(:,:,2*(i-1)+1) = zeros(n,n);
    
    cx(:,:,i) = repmat([1:n]',1,n);
    cy(:,:,i) = repmat([1:n],n,1);
    cz(:,:,i) = (2*(i-1)+1)*ones(n,n);
    
    cg(:,:,i) = v(:,:,2*(i-1)+1);
    
end

cx = (2*cx-(max(cx(:))+min(cx(:))))/( max(cx(:))-min(cx(:)) );
cy = (2*cy-(max(cy(:))+min(cy(:))))/( max(cy(:))-min(cy(:)) );
cz = (2*cz-(max(cz(:))+min(cz(:))))/( max(cz(:))-min(cz(:)) );



%%%%% rearrange observed values and construct mask
xd = cx;
yd = cy;
zd = cz;
pd = cg;


xd = xd(:); yd = yd(:); zd = zd(:); pd = pd(:); pd(find(isnan(pd))) = 0;

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
maskvec = mask(:);

y = observed;

figure, view3Dmodified(v0)
figure, view3Dmodified(y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu1 = 50;
mu2 = 50;
alpha = 100;
opts = struct('maxiters',200,'chambolleit',5,'inner_iters',5,'stopcriterion',2,'tol', 5e-2,'verbose',0);
[clean, obj, times, mses] = rayleighReconstruct3D_v1(y, mask, alpha, mu1, mu2, opts);
clean = reshape(clean,n,n,n);

clean = clean/max(clean(:));
mse = norm(clean(:)-v0(:),2)^2/numel(clean);
mae = sum(abs(clean(:)-v0(:)))/numel(v0);

fprintf('Proposed: Iters: %g, CPU time: %g seconds, MAE: %g\n', length(times)-1,times(end),mae)


figure, view3Dmodified(clean)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alfa = 2.5; % regularization prior
nIter = 100;
t0 = cputime;
[mlr,a,b,xn,yn,zn] = maxLh_interp(cx, cy, cz, cg, n);
[sig, iters_rltv] = usDenoising_v3(mlr,alfa,nIter, 5e-5);
sig = sig./max(sig(:));

tfinal = cputime-t0;
mse_rltv = norm(sig(:)-v0(:),2)^2/numel(sig);
mae_rltv = sum(abs(sig(:)-v0(:)))/numel(v0);

fprintf('RLTV 3D: %d iters, %g seconds, MAE; %g\n', iters_rltv, tfinal, mae_rltv)

figure, view3Dmodified(sig)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

