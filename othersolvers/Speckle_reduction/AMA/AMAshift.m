% Speckle reducction for SAR(synthetic aperture radar) system
% [Related Paper]
% H. Woo and S. Yun, "Alternating minimization algorithm for speckle reduction with a shifting technique",
% IEEE Trans. on Image Processing, 21(4):1701-1714, April 2012.

% Code written by Hyenkyun Woo (hyenkyun@gmail.com, hyenkyun.woo@cc.gatech.edu)
% Version : 2012-10-13-ver1.0


function [uout,n,time,itime]=AMAshift(f,opts)

Gtime = tic;

% f is shifted data...

[M,N]=size(f);

% Initial condition
if opts.model == 1 || opts.model == 2  % I-div model
maxb = opts.max+opts.shift;
minb = opts.shift;
u0 = min(max(f,minb),maxb); 
un = f; 
else % exp model
maxb = log(opts.max+opts.shift);
minb = log(opts.shift);
un = min(max(f,opts.shift),opts.max+opts.shift); 
f = log(f);
u0 = min(max(f,minb),maxb);
end    

n=0; 
p = zeros(M+1,N+1,2);  
z = zeros(M+1,N+1,2);
condition=1;
itime = 0;
VecGeneralParameters = [ M; N; opts.alpha; opts.lambda; minb; maxb; opts.model;];
%% Main Loop
while (condition)
              
        n=n+1; 
    
         Tstart = tic; %%%% Time check    

         if(opts.model==2 || opts.model==4) % acc model
          if opts.L == 1  
            MAXn = 150;
          elseif opts.L == 3  
            MAXn = 100;
          end 
          VecGeneralParameters(3) = opts.alpha*10^(0.3*max((MAXn-n)/MAXn,0));  
         end
        
         % update u,z,p by running one iteration AMA
         [u0,z,p] = AMATV_mex(single(u0),single(z),single(p),single(f),single(VecGeneralParameters));
                      
         itime = itime + toc(Tstart); %%%% Time check 

        %stopping condition 
        if(n>1)
        up = un;
        if(opts.model == 3 || opts.model == 4)
         un = exp(u0);   
        else      
         un = u0;   
        end
        relmse = norm(up - un,'fro')/norm(up,'fro');
        condition=(n<opts.nOuter && relmse>opts.xTol);
        end
        
end
uout = un-opts.shift; 
time = toc(Gtime);