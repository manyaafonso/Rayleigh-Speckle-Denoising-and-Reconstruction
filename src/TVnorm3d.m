function y = TVnorm3d(x)

[M, N, L] = size(x);

tmp = zeros(M,N,L);
tmp(:,1:N-1,:) = x(:,2:N,:);
tmp(:,N,:) = x(:,N,:);

y1 = tmp - x;

tmp = zeros(M,N,L);
tmp(1:M-1,:,:) = x(2:M,:,:);
tmp(M,:,:) = x(M,:,:);

y2 = tmp - x;

tmp = zeros(M,N,L);
tmp(:,:,1:L-1) = x(:,:,2:L);
tmp(:,:,L) = x(:,:,L);

y3 = tmp - x;

% diff = zeros(M,N,L,3);
% diff(:,:,:,1) = y1;
% diff(:,:,:,2) = y2;
% diff(:,:,:,3) = y3;

y = sqrt( sum( y1(:).^2 + y2(:).^2 + y3(:).^2 ) );
