function [cube,np, dimRange] = distributionStep(xd, yd, zd, pd, n, lnn, grid)
% interpolates 3D observations onto a regularized maximum
% likelihood (maxLh) volume estimate (rv)
%
% inputs:   xd, yd, zd - coordinates of each plane
%           pd - observations in each plane
%           n - no. elements of the new volume in each direction (resolution)
%           lnn - Linear neighbourhood number of elements that forms the support region
%                   that influences a voxel. If 1, we get 4 surrounding
%                   cubes. If 2, we get 64... etc
%           grid - string that controls the grid shape, cubic [x,x,x] or
%                   paralelogram [x,x,2x]
%
% outputs:  cube - Array of size [n^3, 4, k] storing the positions and intensity
%                   of each pixel in each voxel support region.
%                   n^3 -> number of final voxels;
%                   4 -> size of vector [x,y,z,pi] storing position (with
%                           axis in the voxel center) and intensity of each
%                           observed pixel;
%                   k -> Maximum number of pixels that exist in the support
%                           region of every voxel
%
%           np - vector where each node voxel equals the Number of neighbour
%               Pixel observations inside the suppport region of each voxel
%
%           dimRange - dimension scale factor

%% -------- Input Management --------------
if ~exist('n','var') || isempty(n), n = 64; end
if ~exist('lnn','var') || isempty(lnn), lnn = 1; end
if ~exist('grid','var') || isempty(grid), grid = 'paralelogram'; end
%% ----------------------------------------
[mx,nx]=size(xd);
[my,ny]=size(yd);
[mz,nz]=size(zd);

xd = xd(:); yd = yd(:); zd = zd(:); pd = pd(:); % 'vectorize' input the matrices
pd(isnan(pd)) = 0;    %Remove error intensity entries

%% the coordinates of each plane are re-scaled to [-1,1] and the scalling
% factors stored for posterior real space plotting.
dimRange =[ (max(xd)-min(xd)), (max(yd)-min(yd)), (max(zd)-min(zd)) ]; %x,y,z dimension scale factors

xn   = (2*xd - min(xd) - max(xd))./dimRange(1);
yn   = (2*yd - min(yd) - max(yd))./dimRange(2);
if strcmp(grid,'cube'),
    zn   = (2*zd - min(zd) - max(zd))./dimRange(3);
else
    zn   = (zd - min(zd) - max(zd))./dimRange(3);
end


%% Setting suport region nodes afected by each pixel observation
% find closest nodes for each observed pixel after obs.
t = ((xn - min(xn))*(n-1))./(max(xn)-min(xn))+1 ;
i_b  = floor(t); %position onf the element bellow
i_a  = (ceil(t)==t & t<n).*(t+1) + (ceil(t)~=t).*ceil(t) + (ceil(t)==t & t>=n).*(t-1); %above
pRelPos(:,1) = t - i_b; %Pixel relative position to the lowest-position neighbour voxel

t = ((yn - min(yn))*(n-1))./(max(yn)-min(yn))+1 ;
j_b = floor(t);
j_a  = (ceil(t)==t & t<n).*(t+1) + (ceil(t)~=t).*ceil(t) + (ceil(t)==t & t>=n).*(t-1);
pRelPos(:,2) = t - j_b;


t = ((zn - min(zn))*(n-1))./(max(zn)-min(zn))+1 ;
k_a  = (ceil(t)==t & t<n).*(t+1) + (ceil(t)~=t).*ceil(t) + (ceil(t)==t & t>=n).*(t-1);
k_b  = floor( t );
pRelPos(:,3) = t - k_b;





%% Initializing matrices
k =5;  % 5->random maximum number of pixel obersvations per voxel
if strcmp(grid,'cube')
    cube = zeros(n^3, 4, k);  % 4->[x y d p] 
    np = zeros(n^3,1);
else
    cube = zeros(2*n^3, 4, k);
    np = zeros(2*n^3,1);
end

for i=1:size(pd,1), % sweeping through all observations
    
    %Location of the 8 surrounding nodes afected by each observation
    ngb = [i_b(i), i_b(i), i_b(i), i_b(i), i_a(i), i_a(i), i_a(i), i_a(i);
           j_b(i), j_b(i), j_a(i), j_a(i), j_b(i), j_b(i), j_a(i), j_a(i);
           k_b(i), k_a(i), k_b(i), k_a(i), k_b(i), k_a(i), k_b(i), k_a(i)]';
   
    for j=1:8, % sweeping through each neighborhood (8 nodes per obs.)
        lindex = sub2ind([n,n,n],ngb(j,1),ngb(j,2),ngb(j,3)); %Linear matrix index
        np(lindex) = np(lindex)+1; 
        cube(lindex,:,np(lindex)) = pRelPos(:,1), pd(i) ];
        
        
%         a(ngb(j,1),ngb(j,2),ngb(j,3))=a(ngb(j,1),ngb(j,2),ngb(j,3))+pd(i)^2/2;
%         b(ngb(j,1),ngb(j,2),ngb(j,3))=b(ngb(j,1),ngb(j,2),ngb(j,3))+1;
%         % b(i,j,k)=0 means that the node (i,j,k) is not covered by any observation
    end
end

% maximum likelihood estimate of regularized grid
mlr = ( (a-eps).*(b~=0)+eps ) ./ ( (b-1).*(b~=0)+1 );


%Plots if no output is selected
if nargout==0,
    figure('name','Distribution Step');
    subplot(1,3,1); mesh( reshape(xd,mx,nx), reshape(yd,my,ny), reshape(zd,mz,nz) ); axis equal
    subplot(1,3,2); mesh( reshape(xn,mx,nx), reshape(yn,my,ny), reshape(zn,mz,nz) ); axis equal
    subplot(1,3,3); mesh( 0.5*dimRange(1)*reshape(xn,mx,nx), 0.5*dimRange(2)*reshape(yn,my,ny), 0.5*dimRange(3)*reshape(zn,mz,nz) );
    axis equal; title('Rescaled');
end
