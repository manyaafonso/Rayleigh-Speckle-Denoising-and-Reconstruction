% --- reconstruction demo --- %
clc; close all; clear
addpath('utils');

% Create the 'original' ultrasound data volume
n = 60;
global v;
v = createObject([n,n,n],10,100,'tube',0);


% Simulate acquisition of 2D images
noPlanes = 30; % user-specified number of planes to acquire
planes = arbitrarySlices(noPlanes,v,n,'random');

% create matrices with coordinates and grey-level intensities from
% structure 'planes'
for idx = 1:noPlanes,
    cx(:,:,idx) = planes(idx).x;
    cy(:,:,idx) = planes(idx).y;
    cz(:,:,idx) = planes(idx).z;
    cg(:,:,idx) = planes(idx).g;
end


% %Perform reconstruction onto a regularized volume
% alfa = 1; % regularization prior
% nIter = 100;
% [mlr,a,b,xn,yn,zn] = maxLh_interp(cx, cy, cz, cg, n);
% sig = usDenoising(mlr,alfa,nIter);
% 
% sig = sig./max(sig(:));
% 
% close
% save temp.mat sig n;
% plotVol();


