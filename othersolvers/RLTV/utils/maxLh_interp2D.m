function [mlr,a,b,xn,yn] = maxLh_interp2D(xd, yd, pd, n)
% interpolates 3D observations onto a regularized maximum
% likelihood (maxLh) volume estimate (rv)
%
% inputs:   xd, yd, zd - coordinates of each plane
%           pd - observations in each plane
%           n - no. elements of the new volume in each direction (resolution)
%
% outputs:  mlr - maximum likelihood estimate on a regularized grid (volume)
%           a - matrix where each node voxel is the sum of its squared neighbour
%               pixel observations
%           b - matrix where each node voxel equals the number of neighbour
%               pixel observations
%           xn - x coordinates in [-1,1]
%           yn - y coordinates in [-1,1]
%           zn - y coordinates in [-1,1]



% 'vectorize' the matrices
xd = xd(:); yd = yd(:); pd = pd(:); pd(find(isnan(pd))) = 0;

% the coordinates of each plane are re-scaled to [-1,1] according to the 'rv' dimensions 
xn   = (2*xd - min(xd) - max(xd))./(max(xd)-min(xd));
yn   = (2*yd - min(yd) - max(yd))./(max(yd)-min(yd));

% find closest nodes for each observed pixel after obs.
t = ((xn - min(xn))*(n-1))./(max(xn)-min(xn))+1 ; i_a  = (ceil(t)==t & t<n).*(t+1) + (ceil(t)~=t).*ceil(t) + (ceil(t)==t & t>=n).*(t-1);
t = ((yn - min(yn))*(n-1))./(max(yn)-min(yn))+1 ; j_a  = (ceil(t)==t & t<n).*(t+1) + (ceil(t)~=t).*ceil(t) + (ceil(t)==t & t>=n).*(t-1);

% before obs.
i_b  = floor( ((xn - min(xn))*(n-1))./(max(xn)-min(xn))+1 );
j_b  = floor( ((yn - min(yn))*(n-1))./(max(yn)-min(yn))+1 );

% matrix a - each node is the sum of its squared neighbors
% matrix b - each node equals the number of neighbor observations
a  = zeros(n,n); 
b  = a;
eps = 0;

for i=1:size(pd,1), % sweeping through all observations
    
    %Location of the 8 surrounding nodes afected by each observation
    ngb = [i_b(i),i_b(i),i_a(i),i_a(i);
            j_b(i),j_a(i),j_b(i),j_a(i)]';
   
    for j=1:4, % sweeping through each neighborhood (8 nodes per obs.)

        a(ngb(j,1),ngb(j,2))=a(ngb(j,1),ngb(j,2))+pd(i)^2/2;
        b(ngb(j,1),ngb(j,2))=b(ngb(j,1),ngb(j,2))+1;
        % b(i,j,k)=0 means that the node (i,j,k) is not covered by any observation
    end
end

% maximum likelihood estimate of regularized grid
mlr = ( (a-eps).*(b~=0)+eps ) ./ ( (b-1).*(b~=0)+1 );
