function y = gradient_3d(x)

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

y = zeros(M,N,L,3);
y(:,:,:,1) = y1;
y(:,:,:,2) = y2;
y(:,:,:,3) = y3;

