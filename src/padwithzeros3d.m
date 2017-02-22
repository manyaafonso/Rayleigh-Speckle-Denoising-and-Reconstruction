function y = padwithzeros3d(L,M,N,x,q)
y = zeros(L,M,N);
y(:,:,q) = x;