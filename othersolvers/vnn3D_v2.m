function [estimated, processing_time] = vnn3D_v2(observed, mask, tx, ty, tz, noPlanes)

t0 = cputime;

[M,N,L] = size(observed);


tx = reshape(tx,M,N,noPlanes);
ty = reshape(ty,M,N,noPlanes);
tz = reshape(tz,M,N,noPlanes);


maskvec = mask(:);

y = observed(:);


%%%% voxel nearest neighbour

estimated = reshape(y,M,N,L);

ll = [1; 1; 1];
ul = [M; N; L];

for i = 1:M
    for j = 1:N
        for k = 1:L
            
            if ~mask(i,j,k)
                Pin = [i;j;k];
                
                dist = [];
                proj_pts = [];

                for slice = 1:noPlanes

                    if maskvec(M*N*(Pin(3)-1)+N*(Pin(2)-1)+Pin(1))==0

                    c = [tx(1,1,slice);ty(1,1,slice);tz(1,1,slice)];
                    b = [tx(1,N,slice);ty(1,N,slice);tz(1,N,slice)];
                    a = [tx(M,1,slice);ty(M,1,slice);tz(M,1,slice)];

%                     ll = [tx(1,1,slice);ty(1,1,slice);tz(1,1,slice)];
%                     ul = [tx(n,n,slice);ty(n,n,slice);tz(n,n,slice)];

                    [po, d] = projectPointonPlane(Pin,a,b,c,ll,ul);
                    dist = [dist, d];
                    proj_pts = [proj_pts, po];

                    end

                end

                ind=find(abs(dist)==min(abs(dist)),1,'first');

                vox_value = y( round( M*N*(proj_pts(3,ind)-1)+N*(proj_pts(2,ind)-1)+proj_pts(1,ind) ) );
                estimated(i,j,k) = vox_value;
    
            end
            
        end
    end
end
    
processing_time = cputime-t0;
