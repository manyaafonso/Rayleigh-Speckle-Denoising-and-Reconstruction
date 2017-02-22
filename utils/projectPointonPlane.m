function [po,d] = projectPointonPlane(p,a,b,c,lowerlim,upperlim)

% if length(b) == 3
    n = cross(b-a,c-a);
% else
%     n = cross([b-a; 0],[c-a; 0]);
%     n = n(1:2);
% end
n = n/norm(n);

d = (p-a)'*n;

po = p-d*n;

po = min(po, upperlim);
po = max(po, lowerlim);
