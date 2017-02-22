function [clean, obj, times, f1, f2, mses] = rayleighDenoise3D_v1(varargin)%(y, alpha, mu, x_true)


if nargin < 3
    error('Must specify observed volume y, regularization parameter alpha, and weight mu!')
else
    y = varargin{1};
    alpha = varargin{2};
    mu = varargin{3};
    
    % default parameters
    compute_mse = 0;
    maxiters = 8;
    chambolleit = 5;
    inner_iters = 5;
    verbose = 0;
    stopcriterion = 0;

    if nargin == 4
        %%% if options are provided, use them, else use default parameters.
        opts = varargin{4};
        if isfield(opts,'maxiters')
            maxiters = opts.maxiters;
        end
        if isfield(opts,'chambolleit')
            chambolleit = opts.chambolleit;
        end
        if isfield(opts,'inner_iters')
            inner_iters = opts.inner_iters;
        end
        if isfield(opts,'x_true')
            compute_mse = 1;
            x_true = opts.x_true;
        end
        if isfield(opts,'verbose')
            verbose = opts.verbose;
        end
        if isfield(opts,'stopcriterion')
            stopcriterion = opts.stopcriterion;
            if isfield(opts,'tol')
                tol = opts.tol;
            else
                if stopcriterion
                    error('Must specify tolerance for stopping criterion.');
                end
            end
        end
        

    end
end

[L,M,N] = size(y);

u = 0.5*y.^2;%zeros(L,M,N);
v = u;
d = u;

ysq = y.^2;

f1(1) = TVnorm3d(u);
f2(1) = sum(0.5*ysq(:).*exp(-u(:))+u(:));
obj(1) = alpha*TVnorm3d(u) + sum(0.5*ysq(:).*exp(-u(:))+u(:));
times(1) = 0;
t0 = cputime;
mses = [];

clean_prev = exp(u);

for t = 1:maxiters
    if verbose
        t
    end
        
    r = v+d;
%     for i = 1:L
%         for j = 1:M
%             for k = 1:N
%                 %u(i,j,k) = denoise1(y(i,j,k),mu,v(i,j,k)+v(i,j,k));
%                 u(i,j,k) = denoise1(y(i,j,k),mu,r(i,j,k));
%             end
%         end
%     end
    u = denoise_substep(y, mu, r, inner_iters);
    clean = exp(u);
    
    v = chambolledenoise3d_v4(u-d, alpha/mu, chambolleit );
    %v = projk3d_v2(u-d, alpha/mu, 4 );
    d = d - (u-v);
    
    if compute_mse
        mses(t) = norm(clean(:)-x_true(:),2)^2/numel(x_true);
    end
    
    f1(t+1) = TVnorm3d(u);
    f2(t+1) = sum(0.5*ysq(:).*exp(-u(:))+u(:));
    obj(t+1) = alpha*TVnorm3d(u) + sum(0.5*ysq(:).*exp(-u(:))+u(:));
    times(t+1) = cputime - t0;
    
    if stopcriterion
        
        switch stopcriterion
            case 1
                criterion = abs( (obj(t+1)-obj(t))/obj(t) );
            case 2
                criterion = norm(clean(:)-clean_prev(:))/norm(clean_prev(:));
                clean_prev = clean;
            case 3
                criterion = obj(t+1);
            otherwise
                error('Invalid stopping criterion!')
        end
        if criterion < tol
            if verbose
                fprintf('Convergence reached.\n')
            end
            break;
        end
        
    end
    
    
    
end
