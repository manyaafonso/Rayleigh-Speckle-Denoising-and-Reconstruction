% Speckle reducction for SAR(synthetic aperture radar) system
% [Related Paper]
% H. Woo and S. Yun, "Alternating minimization algorithm for speckle reduction with a shifting technique",
% IEEE Trans. on Image Processing, 21(4):1701-1714, April 2012.

% Code written by Hyenkyun Woo (hyenkyun@gmail.com, hyenkyun.woo@cc.gatech.edu)
% Version : 2012-10-13-ver1.0


clear all
%close all

%% Setup
Para.L     =  [1 3];   % Noise level (SAR System Parameter)
Para.Model =  [1 2 3 4]; %I-div, I-div(adaptive), Exp, Exp(adaptive)
Iorg = imread('lena.png');
Iorg = double(Iorg);  

[M,N]=size(Iorg); %M:y, N:x
opts.Iorg = Iorg;
opts.nOuter = 400; % max iteration
opts.xTol = 3*10^-4; % relative error 
opts.max   = 255;  % Image range
opts.min   = 0;    % image range
opts.shift = 30;    % (T=30) shift given data

if(isunix==1) 
    disp('fast log is used in Linux');
else
    disp('fast log is not used in Windows');
end   
%% AMA with a shifting technique for speckle reduction

for GL = Para.L       
    
% Add gamma noise    
 opts.L = GL;
 
s = zeros(size(Iorg));
for k = 1:opts.L
    s = s + abs(randn(size(Iorg)) + 1i * randn(size(Iorg))).^2 / 2;
end
Bn = Iorg.*(s/opts.L);    

Bne = Bn + opts.shift;    

   % I-div model
   opts.model = 1;
        if(opts.L==1)
         opts.lambda = 1;             opts.alpha = 15.0*10^-5;  
         [L1out_idiv.Img,L1out_idiv.iter,L1out_idiv.time,L1out_idiv.itime]=AMAshift(Bne,opts);
        else %L=3
         opts.lambda = 1.3/opts.L;    opts.alpha = 24.0*10^-5; 
         [L3out_idiv.Img,L3out_idiv.iter,L3out_idiv.time,L3out_idiv.itime]=AMAshift(Bne,opts);
        end   
     

   % I-div model (adaptive)
   opts.model = 2;
        if(opts.L==1)
         opts.lambda = 1;             opts.alpha = 15.0*10^-5;   
         [L1out_idiva.Img,L1out_idiva.iter,L1out_idiva.time,L1out_idiva.itime]=AMAshift(Bne,opts);
        else %L=3
         opts.lambda = 1.3/opts.L;    opts.alpha = 24.0*10^-5;   
         [L3out_idiva.Img,L3out_idiva.iter,L3out_idiva.time,L3out_idiva.itime]=AMAshift(Bne,opts);
        end 

        
   %exponential model 
   opts.model = 3;
        if(opts.L==1)
         opts.lambda = 1;          opts.alpha = 43*10^-3;  
         [L1out_exp.Img,L1out_exp.iter,L1out_exp.time,L1out_exp.itime]=AMAshift(Bne,opts);
        else %L=3
         opts.lambda = 1.3/opts.L;  opts.alpha = 60*10^-3;  
         [L3out_exp.Img,L3out_exp.iter,L3out_exp.time,L3out_exp.itime]=AMAshift(Bne,opts);
        end   
  
        
   %exponenital model (adaptive)
   opts.model = 4;
        if(opts.L==1)
         opts.lambda = 1;            opts.alpha = 43*10^-3; 
         [L1out_expa.Img,L1out_expa.iter,L1out_expa.time,L1out_expa.itime]=AMAshift(Bne,opts);
        else %L=3
         opts.lambda = 1.3/opts.L;   opts.alpha = 60*10^-3;   
         [L3out_expa.Img,L3out_expa.iter,L3out_expa.time,L3out_expa.itime]=AMAshift(Bne,opts);
        end   
    

end % L


%% plot
figure('name','Alternating minimization algorithm for speckle reduction with a shifting technique(IEEE TIP 2012)'); 

subplot(261); 
imshow(Iorg,[0 255]), title(' (a) Org ');  
subplot(262); 
imshow(Bn,[0 255]), title(' (b) Noisy(L=1) ');

subplot(263); 
MSE=norm(abs(Iorg-L1out_idiv.Img),'fro')^2/(M*N);
psnr=10.*log10(255^2/MSE);
imshow(L1out_idiv.Img,[0 255]), title(['(c) I-DIV, ',num2str(psnr,'%2.2f'),'dB, ',num2str(L1out_idiv.itime,'%2.2f'),'s,',num2str(L1out_idiv.iter,'%2d')]);    

subplot(264);   
MSE=norm(abs(Iorg-L1out_idiva.Img),'fro')^2/(M*N);
psnr=10.*log10(255^2/MSE);
imshow(L1out_idiva.Img,[0 255]), title(['(d) I-DIV acc, ',num2str(psnr,'%2.2f'),'dB, ',num2str(L1out_idiva.itime,'%2.2f'),'s,',num2str(L1out_idiva.iter,'%2d')]);    

subplot(265); 
MSE=norm(abs(Iorg-L1out_exp.Img),'fro')^2/(M*N);
psnr=10.*log10(255^2/MSE);
imshow(L1out_exp.Img,[0 255]), title(['(e) EXP, ',num2str(psnr,'%2.2f'),'dB, ',num2str(L1out_exp.itime,'%2.2f'),'s,',num2str(L1out_exp.iter,'%2d')]);    

subplot(266);   
MSE=norm(abs(Iorg-L1out_expa.Img),'fro')^2/(M*N);
psnr=10.*log10(255^2/MSE);
imshow(L1out_expa.Img,[0 255]), title(['(f) EXP acc, ',num2str(psnr,'%2.2f'),'dB, ',num2str(L1out_expa.itime,'%2.2f'),'s,',num2str(L1out_expa.iter,'%2d')]);    


subplot(268); 
imshow(Bn,[0 255]), title(' (g) Noisy(L=3) ');

subplot(269); 
MSE=norm(abs(Iorg-L3out_idiv.Img),'fro')^2/(M*N);
psnr=10.*log10(255^2/MSE);
imshow(L3out_idiv.Img,[0 255]), title(['(f) I-DIV, ',num2str(psnr,'%2.2f'),'dB, ',num2str(L3out_idiv.itime,'%2.2f'),'s,',num2str(L3out_idiv.iter,'%2d')]);    

subplot(2,6,10);   
MSE=norm(abs(Iorg-L3out_idiva.Img),'fro')^2/(M*N);
psnr=10.*log10(255^2/MSE);
imshow(L3out_idiva.Img,[0 255]), title(['(g) I-DIV acc, ',num2str(psnr,'%2.2f'),'dB, ',num2str(L3out_idiva.itime,'%2.2f'),'s,',num2str(L3out_idiva.iter,'%2d')]);    

subplot(2,6,11); 
MSE=norm(abs(Iorg-L3out_exp.Img),'fro')^2/(M*N);
psnr=10.*log10(255^2/MSE);
imshow(L3out_exp.Img,[0 255]), title(['(f) EXP, ',num2str(psnr,'%2.2f'),'dB, ',num2str(L3out_exp.itime,'%2.2f'),'s,',num2str(L3out_exp.iter,'%2d')]);    

subplot(2,6,12);   
MSE=norm(abs(Iorg-L3out_expa.Img),'fro')^2/(M*N);
psnr=10.*log10(255^2/MSE);
imshow(L3out_expa.Img,[0 255]), title(['(g) EXP acc, ',num2str(psnr,'%2.2f'),'dB, ',num2str(L3out_expa.itime,'%2.2f'),'s,',num2str(L3out_expa.iter,'%2d')]);    
