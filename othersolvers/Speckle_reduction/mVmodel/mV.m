% Speckle reducction for SAR(synthetic aperture radar) system
% [Related Paper]
% S. Yun and H. Woo, "A new multiplicative denoising variational model based on m-th root transformation",
% IEEE Trans. on Image Processing, 21(5):2523-2533, May 2012.

% Code written by Hyenkyun Woo (hyenkyun@gmail.com, hyenkyun.woo@cc.gatech.edu)
% Version : 2012-10-13-ver1.0

function [uout,n,time,itime]=mV(f,opts)

Gtime = tic;
% f is shifted data...
[M,N]=size(f);
maxb = (opts.max+opts.shift)^(1/opts.gsq);
minb = opts.shift;

%% Initial condition
u0 = ones(M,N)*maxb;
un = opts.max+opts.shift; 
n=0; 
p = zeros(M+1,N+1,2);  
z = zeros(M+1,N+1,2);
condition=1;
itime = 0;
VecGeneralParameters = [ M; N; opts.alpha; opts.lambda; opts.gsq; minb; maxb; opts.rho;];
F = opts.gsq*f;
%% Main Loop
while (condition)
               
        n=n+1; 
               
        % update u,z,p by running one iteration LPAMA
        Tstart = tic; %%%% Time check   
        [u0,z,p] = PAMATV_mex(single(u0),single(z),single(p),single(F),single(VecGeneralParameters));                
        itime = itime + toc(Tstart); %%%% Time check
 
        %stopping condition 
        up = un;
        un = u0.*u0;
        if(opts.gsq==4)  un = un.*un;  end
        if(opts.gsq==8)  un = (un.*un); un = (un.*un); end
        if(opts.gsq==32) un = (un.*un); un = (un.*un); un = (un.*un); un = (un.*un); end            
        relmse = norm(up - un,'fro')/norm(up,'fro');
        condition=(n<opts.nOuter && relmse>opts.xTol); 
        
end
uout = un-opts.shift; 
time = toc(Gtime);