function y = divergence3d(x1, x2, x3)

[M, N, L] = size(x1);

y = zeros(M,N,L);

tmp1 = zeros(M,N,L);
tmp2 = tmp1;
tmp1(:,1:N-1,:) = x1(:,1:N-1,:);
tmp2(:,2:N,:) = x1(:,1:N-1,:);

y = y + tmp1-tmp2;

tmp1 = zeros(M,N,L);
tmp2 = tmp1;
tmp1(1:M-1,:,:) = x2(1:M-1,:,:);
tmp2(2:M,:,:) = x2(1:M-1,:,:);

y = y + tmp1-tmp2;

tmp1 = zeros(M,N,L);
tmp2 = tmp1;
tmp1(:,:,1:L-1) = x3(:,:,1:L-1);
tmp2(:,:,2:L) = x3(:,:,1:L-1);

y = y + tmp1-tmp2;

