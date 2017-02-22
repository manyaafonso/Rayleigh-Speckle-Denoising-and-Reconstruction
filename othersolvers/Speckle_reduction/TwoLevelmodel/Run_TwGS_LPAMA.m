% Speckle reducction for SAR(synthetic aperture radar) system
% [Related Paper]
% M. Kang, S. Yun, and H. Woo, "Two-level convex relaxed variational model for multiplicative denoising", 
% submitted to SIAM Journal of Imagining Sciences


% Code written by Hyenkyun Woo (hyenkyun@gmail.com, hyenkyun.woo@cc.gatech.edu)
% Version : 2012-10-13-ver1.0

% Readme:
% Run_TwGS_LPAMA.m : main execution program file
% ...mex : main lib file

clear all
close all

%% Setup
Para.L     =  [1 3];   % Noise level (SAR System Parameter)
Para.mSQ   =  [2 4];     % m-th root transformed variational model

Iorg = imread('lena.png');
Iorg = double(Iorg);  

[M,N]=size(Iorg); %M:y, N:x

opts.nOuter = 500; % max iteration
opts.xTol = 10^-3; % relative error 
opts.dual = 4;     % Concave Dual
opts.KMAX = 5;     % number of inner iteration for dual 
opts.max   = 255;  % Image range
opts.min   = 0;    % image range
opts.init  = 1;    % init (1: max, 2:mean, 3:proj, 4:min);
opts.shift = 1;    %shift given data

%% Two-level relaxed variational model for multiplicative denoising
for GL = Para.L       
    
% Add gamma noise    
 opts.L = GL;
 
s = zeros(size(Iorg));
for k = 1:opts.L
    s = s + abs(randn(size(Iorg)) + 1i * randn(size(Iorg))).^2 / 2;
end
Bn = Iorg.*(s/opts.L);    

Bne = Bn + opts.shift;    
   
for GmSQ   = Para.mSQ
    
% Set m-th root
 opts.gsq = GmSQ; 

 switch GmSQ
 
     case 2 % 2-th root
         
        if(opts.L==1)
         opts.lambda = 0.2;          opts.alpha = 0.03;   opts.rho = 10;
        else %L=3
         opts.lambda = 0.28/opts.L;  opts.alpha = 0.03;   opts.rho = 10; 
        end
        
        [out_sq2.Img,out_sq2.iter,out_sq2.time,out_sq2.itime]=TwLmV(Bne,opts);

     case 4 % 4-th root
         
        if(opts.L==1)
         opts.lambda = 1.3;         opts.alpha = 1.0;    opts.rho = 0.3;   
        else %L=3
         opts.lambda = 2/opts.L;    opts.alpha = 1.0;    opts.rho = 0.3; 
        end 
    
        [out_sq4.Img,out_sq4.iter,out_sq4.time,out_sq4.itime]=TwLmV(Bne,opts);
        
    
 end
 
end % mSQ  
   if GL==1
     L1out_sq2 = out_sq2;
     L1out_sq4 = out_sq4;
   elseif GL==3
     L3out_sq2 = out_sq2;
     L3out_sq4 = out_sq4;
   else
     disp('err');  
   end    
end % L


%% plot
figure('name','Two-level convex relaxed variational model for multiplicative denoising(submitted to SIAM JIS)'); 

subplot(241); 
imshow(Iorg,[0 255]), title(' (a) Org ');  
subplot(242); 
imshow(Bn,[0 255]), title(' (b) Noisy(L=1) ');

subplot(243); 
MSE=norm(abs(Iorg-L1out_sq2.Img),'fro')^2/(M*N);
psnr=10.*log10(255^2/MSE);
imshow(L1out_sq2.Img,[0 255]), title(['(c) TwL-2V, ',num2str(psnr,'%2.2f'),'dB, ',num2str(L1out_sq2.itime,'%2.2f'),'s,',num2str(L1out_sq2.iter,'%2d')]);    

subplot(244);   
MSE=norm(abs(Iorg-L1out_sq4.Img),'fro')^2/(M*N);
psnr=10.*log10(255^2/MSE);
imshow(L1out_sq4.Img,[0 255]), title(['(d) TwL-4V, ',num2str(psnr,'%2.2f'),'dB, ',num2str(L1out_sq4.itime,'%2.2f'),'s,',num2str(L1out_sq4.iter,'%2d')]);    

subplot(246); 
imshow(Bn,[0 255]), title(' (e) Noisy(L=3) ');

subplot(247); 
MSE=norm(abs(Iorg-L3out_sq2.Img),'fro')^2/(M*N);
psnr=10.*log10(255^2/MSE);
imshow(L3out_sq2.Img,[0 255]), title(['(f) TwL-2V, ',num2str(psnr,'%2.2f'),'dB, ',num2str(L3out_sq2.itime,'%2.2f'),'s,',num2str(L3out_sq2.iter,'%2d')]);    

subplot(248);   
MSE=norm(abs(Iorg-L3out_sq4.Img),'fro')^2/(M*N);
psnr=10.*log10(255^2/MSE);
imshow(L3out_sq4.Img,[0 255]), title(['(g) TwL-4V, ',num2str(psnr,'%2.2f'),'dB, ',num2str(L3out_sq4.itime,'%2.2f'),'s,',num2str(L3out_sq4.iter,'%2d')]);    
