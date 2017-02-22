function y = padwithzeros(x,mask)
y = zeros(size(mask));
y(mask==1) = x;