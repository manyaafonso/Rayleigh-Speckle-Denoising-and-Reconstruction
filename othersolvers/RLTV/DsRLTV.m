function [sig, eta, Ed, Ep, t] = DsRLTV(y, a, opt)
%% DsRLTV   - Total Variation De-Speckling using a Rayleigh Observation
% Model and Log Euclidean Priors
% 
% jose seabra ISR-IST Portugal  
% Updated: 21.07.2011    
% ------------------------------------------------------------------------
% Usage:  [sig, eta, Ed, Ep, t] = DsRLTV(your_image, 1, 1)
%
% INPUTS: 
% y        - double format image
% a(alpha) - regularization parameter:
%            1-10: manual {noisier->smoother}
%            -1,-2,-3: auto {alpha Benford, alpha Optimal, alpha Filtered}
% opt      - visualization of results:
%            0: no profiles are plotted
%            1: display image profile, energy and gradient plots        
% ------------------------------------------------------------------------
% OUTPUT: 
% sig      - de-noised (reconstructed) image
% eta      - speckle field estimate
% Ed       - Energy of data term
% Ep       - Energy of prior term
% t        - elapsed time
%% ------------------------------------------------------------------------
clc
close all

% check (and correct) image type
if size(y,3) ~= 1,
    y = rgb2gray(y);  
end
y = double(y);


% ---- Newton method [cooling + predictor] ----
eps_newton = 1;
number_solved_problems = 0;
tic

[n,m] = size(y);
ve = [];
eps = 1e-3;
A = 0.5*y.^2 + eps;
x = log(A);
ev = [];
ng=[];

%% Set the Regularization Hyperparameter (alpha)
if a<0 % automatic alpha
        [gx,gy]=gradient(x);     G=sqrt(gx.^2+gy.^2);
        G(G==0)=eps;
        if      a == -1,   alpha=1./G;                              %Alpha Benford
        elseif  a == -2,   alpha=1/mean(G(:));                      %Alpha Optimal
        elseif  a == -3,   alpha=1./conv2(G, ones(3)./3^2 ,'same'); %Alpha Filtered
        end
        
elseif a == 0, alpha = 0; %the solution is the Maximum Likelihood estimate
else
    alpha = a; % manual alpha
end

%% compute the Energy function (f)
dif_j = diff(x')'; dif_j2 = dif_j.*dif_j;   aux_j = [zeros(n,1), dif_j];
dif_i = diff(x);   dif_i2 = dif_i.*dif_i;   aux_i = [zeros(1,m); dif_i];

aux_j2 = [zeros(n,1), dif_j2]; 
aux_i2 = [zeros(1,m); dif_i2];

aux_ij = aux_i + aux_j;
grad_ij = sqrt(aux_i2 + aux_j2 + eps_newton); 

f = sum(sum(A.*exp(-x) + x + alpha.*grad_ij));

%% compute Gradient of f
grad_j = [grad_ij(:,2:end), grad_ij(:,end)];
grad_i = [grad_ij(2:end,:); grad_ij(end,:)];

dif_j = [-dif_j, zeros(n,1)];
dif_i = [-dif_i; zeros(1,m)];

g = -A.*exp(-x) + 1 + alpha*(aux_ij./grad_ij + (dif_j./grad_j) + (dif_i./grad_i));

% indices (for heptagonal hessian matrix)
indxs_i = [ 1:n*m   2:n*m       1:n*m-1    n:n*m       1:n*m-n+1   n+1:n*m     1:n*(m-1)]';
indxs_j = [ 1:n*m   1:n*m-1     2:n*m      1:n*m-n+1   n:n*m       1:n*(m-1)   n+1:n*m]';


%% optimization loop
iter = 0;
while ((iter < 2e2) && (eps_newton >= 1e-6))
  
    % compute gradient of f
    dif_j = diff(x')'; dif_j2 = dif_j.*dif_j;   aux_j = [zeros(n,1), dif_j];
    dif_i = diff(x);   dif_i2 = dif_i.*dif_i;   aux_i = [zeros(1,m); dif_i];
    aux_j2 = [zeros(n,1), dif_j2]; 
    aux_i2 = [zeros(1,m); dif_i2];
    
    aux_ij = aux_i + aux_j;

    dif1 = -[dif_j, zeros(n,1)];
    dif2 = -[aux_i(:,2:end), zeros(n,1)];
    
    grad_ij = sqrt(aux_i2 + aux_j2 + eps_newton); 
    grad_j = [grad_ij(:,2:end), grad_ij(:,end)];
    grad_i = [grad_ij(2:end,:); grad_ij(end,:)];
    
    dif_j = [-dif_j, zeros(n,1)];
    dif_i = [-dif_i; zeros(1,m)];

    g = -A.*exp(-x) + 1 + alpha*(aux_ij./grad_ij + (dif_j./grad_j) + (dif_i./grad_i));
    
    ng = [ng, norm(g)];
    
    % ---- check stopping criterion ----
    if norm(g(:)) < 1e-4
        
            number_solved_problems = number_solved_problems + 1;                    
            new_eps_newton = eps_newton*0.3;    
   
            % predictor code
            if number_solved_problems >= 2
                estimated_speed = (x - x_old_predictor)/(eps_newton - eps_old);
            
                x_old_predictor = x;
                eps_old = eps_newton;
                
                x = x+(new_eps_newton - eps_newton)*estimated_speed;                                
            else
                x_old_predictor = x;
                eps_old = eps_newton;
            end;
                                                
            eps_newton = new_eps_newton;
       
            % compute f  
            dif_j = diff(x')'; dif_j2 = dif_j.*dif_j;   aux_j = [zeros(n,1), dif_j];
            dif_i = diff(x);   dif_i2 = dif_i.*dif_i;   aux_i = [zeros(1,m); dif_i];

            aux_j2 = [zeros(n,1), dif_j2]; 
            aux_i2 = [zeros(1,m); dif_i2];        
            aux_ij = aux_i + aux_j;

            grad_ij = sqrt(aux_i2 + aux_j2 + eps_newton); 

            dif1 = -[dif_j, zeros(n,1)];
            dif2 = -[aux_i(:,2:end), zeros(n,1)];            

            
                % ---- alpha ----
                if a<0 % automatic alpha
        
                    [gx,gy]=gradient(x);     G=sqrt(gx.^2+gy.^2);        
        
                    G(G==0)=eps;
        
                    if      a == -1,   alpha=1./G;                              %Alpha Benford
                    elseif  a == -2,   alpha=1/mean(G(:));                      %Alpha optimo
                    elseif  a == -3,   alpha=1./conv2(G, ones(3)./3^2 ,'same'); %Alpha Filtrado
                    end
                elseif a == 0, alpha = 0; %the solution is the Maximum Likelihood estimate
    
                else

                    alpha = a; % manual alpha
    
                end
            
            
            f(end) = sum(sum(A.*exp(-x) + x + alpha.*grad_ij));

            % compute gradient of f 
            grad_j = [grad_ij(:,2:end), grad_ij(:,end)];
            grad_i = [grad_ij(2:end,:); grad_ij(end,:)];
    
            dif_j = [-dif_j, zeros(n,1)];
            dif_i = [-dif_i; zeros(1,m)];

            g = -A.*exp(-x) + 1 + alpha*(aux_ij./grad_ij + (dif_j./grad_j) + (dif_i./grad_i));
            
    end
    
    % compute hessian               
    aux_i2 = [aux_i2(:,2:end), zeros(n,1)]; aux_i2(:,1:end-1) = aux_i2(:,1:end-1) + eps_newton; 
    aux_j2 = [aux_j2(2:end,:); zeros(1,m)]; aux_j2(1:end-1,:) = aux_j2(1:end-1,:) + eps_newton;

    aux_ij = 2.*eps_newton.*ones(n,m); aux_ij = aux_ij + aux_j.^2 + aux_i.^2 - 2.*aux_j.*aux_i; 
    aux_ij(1,:) = eps_newton;   aux_ij(:,1) = eps_newton;  aux_ij(1,1) = 0;

    h0 = A.*exp(-x) + alpha*(aux_ij./grad_ij.^3 + aux_i2./grad_j.^3 + aux_j2./grad_i.^3);
    
    a_i = [aux_i(:,2:end), zeros(n,1)];
    a_j = [zeros(1,m); -dif_j(1:end-1,:)]; a_j(:,end) = 0;
    
    d2 = alpha.*((a_i.*(a_j - aux_i) - eps_newton)./grad_j.^3)'; d2=d2(1:end-1,:);
    
    b_j = [ aux_j(2:end,:); zeros(1,m)];
    b_i = [zeros(n,1), -dif_i(:,1:end-1)];

    d1 = alpha.*((b_j.*(b_i - aux_j) - eps_newton)./grad_i.^3)'; d1(:,end) = 0; 
    d3 = -alpha * (dif1.*dif2)./(dif1.^2 + dif2.^2 + eps_newton).^1.5;
        
    [l c]=size(h0);  h0=reshape(h0,1,l*c);
    [l c]=size(d1);  d1=reshape(d1',1,l*c); d1 = d1(1:end-1);
    [l c]=size(d2);  d2=reshape(d2',1,l*c);
    [l c]=size(d3);  d3=reshape(d3,1,l*c); d3 = d3(1:end-(n-1));
    
    H = sparse(indxs_i,indxs_j,[ h0' ; d1' ; d1' ; d3'; d3'; d2' ; d2' ],n*m,n*m);

    % set search direction
    dvec = -H\g(:);    
    d = vec2mat(dvec,n)';
    
    % perform Armijo rule
    s  = 1;
    xs = x + s.*d;
    f0 = f(end);
    n_a = 0;
    
    dif_js = diff(xs')'; dif_js2 = dif_js.*dif_js;   aux_js = [zeros(n,1), dif_js];
    dif_is = diff(xs);   dif_is2 = dif_is.*dif_is;   aux_is = [zeros(1,m); dif_is];

    aux_js2 = [zeros(n,1), dif_js2]; 
    aux_is2 = [zeros(1,m); dif_is2];
    aux_ijs = aux_is + aux_js;

    grad_ijs = sqrt(aux_is2 + aux_js2 + eps_newton);     
    slope = (1e-4).*g(:)'*dvec; 
 
    while sum(sum(A.*exp(-xs) + xs + alpha.*grad_ijs)) > f0+s*slope,
        s = 0.5*s;
        xs = x+s.*d;
        n_a = n_a+1;
        if n_a > 200
            fprintf('Armijo problems...sparse (predictor)'); keyboard;
        end;
                
        dif_js = diff(xs')'; dif_js2 = dif_js.*dif_js;   aux_js = [zeros(n,1), dif_js];
		dif_is = diff(xs);   dif_is2 = dif_is.*dif_is;   aux_is = [zeros(1,m); dif_is];

		aux_js2 = [zeros(n,1), dif_js2]; 
		aux_is2 = [zeros(1,m); dif_is2];

		aux_ijs = aux_is + aux_js;
		grad_ijs = sqrt(aux_is2 + aux_js2 + eps_newton);
        
    end;
	x = xs;
    
    % ---- alpha ----
    if a<0 % automatic alpha
        [gx,gy]=gradient(x);     G=sqrt(gx.^2+gy.^2);        
        G(G==0)=eps;
        if      a == -1,   alpha=1./G;                              %Alpha Benford
        elseif  a == -2,   alpha=1/mean(G(:));                     %Alpha optimo
        elseif  a == -3,   alpha=1./conv2(G, ones(3)./3^2 ,'same'); %Alpha Filtrado
        end
    elseif a == 0, alpha = 0; %the solution is the Maximum Likelihood estimate
    else
        alpha = a; % manual alpha
    end
    
    f = [ f , sum(sum(A.*exp(-x) + x + alpha.*grad_ijs)) ];  
    iter = iter+1; ev = [ev iter]; clc;
    disp(ev);
    
    
    sig= sqrt(exp(x));
    eta = y./(2.*sqrt(sig));
    
    %% plot results at each iteration
% % %     figure(100); mesh(sig); colormap gray; axis off; axis square; drawnow; 
% % %     figure(101); imagesc(log(eta/2+1)); colormap gray; axis off; axis square; drawnow; 
% % %     ve = [ve, log(eta(end/2,end/2)/2+1)];
% % %     figure(102); plot(ve,'x'); drawnow,    
end
% end of loop

t = toc; % check elapsed time

f9 = f;

Ed = sum(sum(A.*exp(-x) + x)); % Energy of Data Term
Ep = sum(sum(grad_ijs));       % Energy of Prior Term

ng9 = ng;

sig= sqrt(exp(x));
eta = y./(2.*sqrt(sig));    % speckle field estimate


%% plot final results
clc
    figure; imagesc(sqrt(A)); colormap gray; axis off; axis square; hold on; colorbar
    title('original image');
    
    figure; imagesc(sig); colormap gray; axis off; axis square; colorbar
    title('de-speckled image');
    
    figure; imagesc(eta); colormap gray; axis off; axis square; colorbar
    title('speckle image');
    
    
if opt==1, 
    
    figure; semilogx(f9,'k-'); grid on; title('Energy function');
    hold on; xlabel('iteration'); hold on; ylabel('E'); axis tight;
    
    figure; semilogx(ng9,'k-'); grid on; title('Gradient norm');
    title(sprintf('Gradient norm - convergence in: %i iterations, runtime: %f s',iter, t));
    hold on; xlabel('iteration'); hold on; ylabel('\nabla(E)'); axis tight;
      
    figure; plot(diag(sig),'LineWidth',2);
    title('Image profiles'); hold on; plot(diag(sqrt(A)),'m'); axis tight;
end