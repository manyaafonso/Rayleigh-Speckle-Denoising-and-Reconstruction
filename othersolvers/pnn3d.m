function [estimated, processing_time] = pnn3d(observed, mask, WinSize)

t0 = cputime;

[M,N,L] = size(observed);

estimated = observed;

for i = 1:M
    for j = 1:N
        for k = 1:L
            if ~mask(i,j,k)
                count = 0;
                sumNeighbors = 0;
                for p = max(i-WinSize,1):min(i+WinSize,M)
                    for q = max(j-WinSize,1):min(j+WinSize,N)
                        for r = max(k-WinSize,1):min(k+WinSize,L)
                            if mask(p,q,r)
                                sumNeighbors = sumNeighbors + observed(p,q,r);
                                count = count + 1;
                            end
                        end
                    end
                end
                if count>0
                    estimated(i,j,k) = sumNeighbors/count;
                end
            end
        end
    end
end

processing_time = cputime-t0;
