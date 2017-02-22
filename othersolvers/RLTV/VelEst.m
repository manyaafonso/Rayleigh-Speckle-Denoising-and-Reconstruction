clc; close all; clear
addpath('.\utils');

%% Create the 'original' ultrasound data volume
n = 128;
global v;
v = createObject(n,10,100,'tube',0);

%% Sample data
noPlanes = 30; % user-specified number of planes to acquire
planes = arbitrarySlices(noPlanes,v,n,'regular');

% create matrices with coordinates and grey-level intensities from
% structure 'planes'
for idx = 1:noPlanes,
    cx(:,:,idx) = planes(idx).x;
    cy(:,:,idx) = planes(idx).y;
    cz(:,:,idx) = planes(idx).z;
    cg(:,:,idx) = planes(idx).g;
end

%% Image reshaping, vertically and horizontally
[m,n,o] = size(cg);
vr_v = reshape(cg,m*n,o);  % Vertical vectorization
for k = 1:o,                % Horizontal vectorization
    aux = cg(:,:,k)';
    vr_h(:,k) = aux(:);
end
clear aux;

%% Perform reconstruction onto a regularized volume
% alfa = 1; % regularization prior
% nIter = 100;
% [mlr,a,b,xn,yn,zn] = maxLh_interp(cx, cy, cz, cg, n);

%% Perform image denoising
% sig = usDenoising(mlr,alfa,nIter);
% sig = sig./max(sig(:)); 

%% Corr evolution on sampled volume
filter_lag = length(filter)-1;
for i = 1:o-filter_lag
    for j =1:filter_lag
        %% vertical Xcorr lags
        [r,lags] = xcorr( vr_v(:,i)-mean(vr_v(:,i)),  vr_v(:,i+j)-mean(vr_v(:,i+j)),'coeff' );
        
        [xcorr_max, ind] = max(r);
        lags_v(i,j) = lags(ind);          %slide to make maximum cross correlation with vertical vectorization
        xcorr_max_v(i,j) = xcorr_max;     %maximum normalized cross correlation
        norm_v(i,j) = norm(vr_v(:,i)-vr_v(:,i+j));
        
        %% Horizontal Xcorr lags
        
        [r,lags] = xcorr( vr_h(:,i)-mean(vr_h(:,i)), vr_h(:,i+j)-mean(vr_h(:,i+j)),'coeff' );
        [xcorr_max, ind] = max(r);
        lags_h(i,j) = lags(ind);
        xcorr_max_h(i,j) = xcorr_max;     %maximum normalized cross correlation
        norm_h(i,j) = norm(vr_h(:,i)-vr_h(:,i+j));
    end
    
    disp(['Frame nº: ', num2str(i), ' of ',num2str(size(vr_v,2))]);
end

%% Show results
%plotVol();

figure('name','test'),
x=1:o-length(filter)-1;
subplot(2,2,1);  plot(x, xcorr_max_h(:,1),'--b',x, xcorr_max_h(:,2),'--r', x, xcorr_max_h(:,3),'--g'); title('Xcorr max Horizontal'); hold on;
subplot(2,2,2);  plot(x, xcorr_max_v(:,1),'--b',x, xcorr_max_v(:,2),'--r', x, xcorr_max_v(:,3),'--g'); title('Xcorr max Vertical'); hold on;
subplot(2,2,3);  plot(x, norm_h(:,1),'--b',x, norm_h(:,2),'r', x, norm_h(:,3),'--g'); title('Norm Horizontal'); hold on;
subplot(2,2,4);  plot(x, norm_v(:,1),'--b',x, norm_v(:,2),'r', x, norm_v(:,3),'--g'); title('Norm Vertical'); hold on;


H = fspecial('gaussian', [200,1], 200);
subplot(2,2,1);  plot(x, imfilter(xcorr_max_h(:,1),H,'replicate'),'k','Linewidth',1.5);
subplot(2,2,2);  plot(x, imfilter(xcorr_max_v(:,1),H,'replicate'),'k','Linewidth',1.5);
subplot(2,2,3);  plot(x, imfilter(norm_h(:,1),H,'replicate'),'k','Linewidth',1.5);
subplot(2,2,4);  plot(x, imfilter(norm_v(:,1),H,'replicate'),'k','Linewidth',1.5);


figure('name','Teste 1'),
subplot(2,2,1); boxplot(xcorr_max_h); title('Xcorr max Horizontal');hold on;
plot(mean(xcorr_max_h,1),'k');
subplot(2,2,2); boxplot(xcorr_max_v); title('Xcorr max Vertical');
subplot(2,2,3); boxplot(norm_h); title('Norm max Horizontal');
subplot(2,2,4); boxplot(norm_v); title('Norm max Vertical');
