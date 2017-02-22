function [estimated, processing_time] = pbmdw(observed, mask, kernel_diameter)

t0 = cputime;

[M,N,L] = size(observed);


WinSize = floor(kernel_diameter/2);

%center = ceil(kernel_diameter/2);%kernel_diameter+1;

kernel = zeros( WinSize+1, WinSize+1, WinSize+1 );

for i = 0:WinSize
    for j = 0:WinSize
        for k = 0:WinSize
            
            if (i+j+k)>0
            %kernel(i,j,k) = exp(-( (i-center)^2+(j-center)^2+(k-center)^2)/2);
                kernel(i+1,j+1,k+1) = 1/( i^2+j^2+k^2 );
            end
        end
    end
end

estimated = observed;

for i = 1:M
    for j = 1:N
        for k = 1:L
            
            if ~mask(i,j,k)
                
                sumWts = 0;
                sumVox = 0;

                for p = max(i-WinSize,1):min(i+WinSize,M)
                    for q = max(j-WinSize,1):min(j+WinSize,N)
                        for r = max(k-WinSize,1):min(k+WinSize,L)
                            
                            if mask(p,q,r)
                                %wt = 1/((p-i)^2 + (q-j)^2 + (r-k)^2);
                                wt = kernel(abs(p-i)+1,abs(q-j)+1,abs(r-k)+1);
                                sumVox = sumVox + wt*observed(p,q,r);
                                sumWts = sumWts + wt;
                                %count = count + 1;
                            end
                            
                        end
                    end
                end
                
                if sumWts>0
                    estimated(i,j,k) = sumVox/sumWts;
                end
                
            end
            
        end
    end
end

processing_time = cputime-t0;
