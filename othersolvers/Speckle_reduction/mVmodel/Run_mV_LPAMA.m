% Speckle reducction for SAR(synthetic aperture radar) system
% [Related Paper]
% S. Yun and H. Woo, "A new multiplicative denoising variational model based on m-th root transformation",
% IEEE Trans. on Image Processing, 21(5):2523-2533, May 2012.

% Code written by Hyenkyun Woo (hyenkyun@gmail.com, hyenkyun.woo@cc.gatech.edu)
% Version : 2012-10-13-ver1.0


clear all
close all

%% Setup
Para.L     =  [1 3];   % Noise level (SAR System Parameter)
Para.mSQ   =  [2 4 8 32];     % m-th root transformed variational model

Iorg = imread('lena.png');

Iorg = double(Iorg);  

[M,N]=size(Iorg); %M:y, N:x

opts.nOuter = 500; % max iteration
opts.xTol = 5*10^-4; % relative error 
opts.max   = 255;  % Image range
opts.min   = 0;    % image range
opts.shift = 1;    %shift given data

%% m-th root transform based variational model for multiplicative denoising
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
         opts.lambda = 0.23;          opts.alpha = 0.1;   opts.rho = 1.0;
        else %L=3
         opts.lambda = 0.3/opts.L;    opts.alpha = 0.1;   opts.rho = 1.4; 
        end   
        [out_sq2.Img,out_sq2.iter,out_sq2.time,out_sq2.itime]=mV(Bne,opts);

     case 4 % 4-th root
        if(opts.L==1)
         opts.lambda = 1.3;         opts.alpha = 2.0;    opts.rho = 0.1;   
        else %L=3
         opts.lambda = 2.0/opts.L;    opts.alpha = 3.5;    opts.rho = 0.1; 
        end 
        [out_sq4.Img,out_sq4.iter,out_sq4.time,out_sq4.itime]=mV(Bne,opts);
        
     case 8 % 8-th root
        if(opts.L==1)
         opts.lambda = 5.0;          opts.alpha = 30.0;   opts.rho = 0.01;
        else %L=3
         opts.lambda = 7.5/opts.L;  opts.alpha = 30.0;   opts.rho = 0.01; 
        end   
        [out_sq8.Img,out_sq8.iter,out_sq8.time,out_sq8.itime]=mV(Bne,opts);
        
     case 32 % 32-th root
        if(opts.L==1)
         opts.lambda = 30.0;          opts.alpha = 600.0;   opts.rho = 0.0004;
        else %L=3
         opts.lambda = 45.0/opts.L;   opts.alpha = 600.0;   opts.rho = 0.0005; 
        end   
        [out_sq32.Img,out_sq32.iter,out_sq32.time,out_sq32.itime]=mV(Bne,opts);     
 end
 
end % mSQ 
   if GL==1
     L1out_sq2  = out_sq2;
     L1out_sq4  = out_sq4;
     L1out_sq8  = out_sq8;
     L1out_sq32 = out_sq32;
   elseif GL==3
     L3out_sq2  = out_sq2;
     L3out_sq4  = out_sq4;
     L3out_sq8  = out_sq8;
     L3out_sq32 = out_sq32;
   else
     disp('err');  
   end   
end % L


%% plot
figure('name','A new muiltiplicative denoising variational model based on mth root transformation(IEEE TIP 2012)'); 

subplot(261); 
imshow(Iorg,[0 255]), title(' (a) Org ');  
subplot(262); 
imshow(Bn,[0 255]), title(' (b) Noisy(L=1) ');

subplot(263); 
MSE=norm(abs(Iorg-L1out_sq2.Img),'fro')^2/(M*N);
psnr=10.*log10(255^2/MSE);
imshow(L1out_sq2.Img,[0 255]), title(['(c) 2-V, ',num2str(psnr,'%2.2f'),'dB, ',num2str(L1out_sq2.itime,'%2.2f'),'s,',num2str(L1out_sq2.iter,'%2d')]);    

subplot(264);   
MSE=norm(abs(Iorg-L1out_sq4.Img),'fro')^2/(M*N);
psnr=10.*log10(255^2/MSE);
imshow(L1out_sq4.Img,[0 255]), title(['(d) 4-V, ',num2str(psnr,'%2.2f'),'dB, ',num2str(L1out_sq4.itime,'%2.2f'),'s,',num2str(L1out_sq4.iter,'%2d')]);    

subplot(265); 
MSE=norm(abs(Iorg-L1out_sq8.Img),'fro')^2/(M*N);
psnr=10.*log10(255^2/MSE);
imshow(L1out_sq8.Img,[0 255]), title(['(e) 8-V, ',num2str(psnr,'%2.2f'),'dB, ',num2str(L1out_sq8.itime,'%2.2f'),'s,',num2str(L1out_sq8.iter,'%2d')]);    

subplot(266);   
MSE=norm(abs(Iorg-L1out_sq32.Img),'fro')^2/(M*N);
psnr=10.*log10(255^2/MSE);
imshow(L1out_sq32.Img,[0 255]), title(['(f) 32-V, ',num2str(psnr,'%2.2f'),'dB, ',num2str(L1out_sq32.itime,'%2.2f'),'s,',num2str(L1out_sq32.iter,'%2d')]);    


subplot(268); 
imshow(Bn,[0 255]), title(' (g) Noisy(L=3) ');

subplot(269); 
MSE=norm(abs(Iorg-L3out_sq2.Img),'fro')^2/(M*N);
psnr=10.*log10(255^2/MSE);
imshow(L3out_sq2.Img,[0 255]), title(['(f) 2-V, ',num2str(psnr,'%2.2f'),'dB, ',num2str(L3out_sq2.itime,'%2.2f'),'s,',num2str(L3out_sq2.iter,'%2d')]);    

subplot(2,6,10);   
MSE=norm(abs(Iorg-L3out_sq4.Img),'fro')^2/(M*N);
psnr=10.*log10(255^2/MSE);
imshow(L3out_sq4.Img,[0 255]), title(['(g) 4-V, ',num2str(psnr,'%2.2f'),'dB, ',num2str(L3out_sq4.itime,'%2.2f'),'s,',num2str(L3out_sq4.iter,'%2d')]);    

subplot(2,6,11); 
MSE=norm(abs(Iorg-L3out_sq8.Img),'fro')^2/(M*N);
psnr=10.*log10(255^2/MSE);
imshow(L3out_sq8.Img,[0 255]), title(['(f) 8-V, ',num2str(psnr,'%2.2f'),'dB, ',num2str(L3out_sq8.itime,'%2.2f'),'s,',num2str(L3out_sq8.iter,'%2d')]);    

subplot(2,6,12);   
MSE=norm(abs(Iorg-L3out_sq32.Img),'fro')^2/(M*N);
psnr=10.*log10(255^2/MSE);
imshow(L3out_sq32.Img,[0 255]), title(['(g) 32-V, ',num2str(psnr,'%2.2f'),'dB, ',num2str(L3out_sq32.itime,'%2.2f'),'s,',num2str(L3out_sq32.iter,'%2d')]);    
