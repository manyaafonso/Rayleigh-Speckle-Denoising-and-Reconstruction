% Paper   : Two-Level convex relaxed variational model for multiplicative denoising
% Authors : Myungjoo Kang, Sangwoon Yun, and Hyenkyun Woo(Corresponding Author)
% Journal : SIAM Journal of Imaging Sciences (submitted)
% Code written by Hyenkyun Woo (hyenkyun@gmail.com, hyenkyun.woo@cc.gatech.edu)
% Version : 2012-10-13-ver1.0

function [uout,n,time,itime]=TwLmV(f,opts)

time = cputime;
% f is shifted data...
[M,N]=size(f);
maxb = (opts.max+opts.shift)^(1/opts.gsq);
minb = opts.shift;

%% Initial condition
switch(opts.init)
    case 1
        u0 = ones(M,N)*maxb; 
        un = ones(M,N)*(opts.max+opts.shift);
    case 2
        meanf = mean(f(:));
        u0 = ones(M,N)*meanf^(1/opts.gsq); 
        un = ones(M,N)*meanf;
    case 3
        u0 = min( max( f.^(1/opts.gsq) , minb) ,maxb); 
        un = min(max(f, minb),maxb);
    case 4
        u0 = ones(M,N)*minb; 
        un = ones(M,N)*opts.shift;
end

n=0; 
ux = u0;
p = zeros(M+1,N+1,2);  
z = zeros(M+1,N+1,2);
condition=1;
itime = 0;
VecGeneralParameters = [ M; N; opts.alpha; opts.lambda; opts.gsq; minb; maxb; opts.rho; opts.dual];
%% Main Loop
while (condition)
        
        n=n+1; 
        
        itimex = cputime; %%%% Time check
        % update u,z,p by running LPAMA 
        [u0,z,p] = PAMATVDual_mex(single(u0),single(ux),single(z),single(p),single(f),single(VecGeneralParameters));          
        % update concave dual variable  
        if(mod(n,opts.KMAX)==0)       
           ux=u0; 
        end    
        itime = itime + (cputime-itimex); %%%% Time check
                
       %% check stopping condition  
        up = un; 
        un = u0.*u0;
        if(opts.gsq==4)
            un = un.*un; 
        end
        relmse = norm(up - un,'fro')/norm(up-ones(M,N),'fro');
        condition=(n<opts.nOuter && relmse>opts.xTol); 
         
end
uout = un-opts.shift; 
time = cputime - time;