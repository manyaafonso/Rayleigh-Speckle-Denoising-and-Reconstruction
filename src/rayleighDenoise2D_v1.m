function [clean, obj, times, f1, f2, mses] = rayleighDenoise2D_v1(varargin)%(y, alpha, mu, x_true)


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
    u = 0.5*y.^2;

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
        if isfield(opts,'init')
            u = opts.init;
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

[M,N] = size(y);

%u = 0.5*y.^2;%zeros(L,M,N);
v = 0*u;
d = 0*u;

ysq = y.^2;

f1(1) = TVnorm(u);
f2(1) = sum(0.5*ysq(:).*exp(-u(:))+u(:));
obj(1) = alpha*TVnorm(u) + sum(0.5*ysq(:).*exp(-u(:))+u(:));
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
    %fprintf('DEBUG 1..\n');
    u = denoise_substep(y, mu, r, inner_iters);
    clean = exp(u);
    %fprintf('DEBUG 2..\n');
    %v = chambolledenoise3d_v2(u-d, alpha/mu, chambolleit );
    v = projk(u-d, alpha/mu, chambolleit );
    %v = projk3d_v2(u-d, alpha/mu, 4 );
    d = d - (u-v);
    %fprintf('DEBUG 3..\n');
    if compute_mse
        mses(t) = norm(clean(:)-x_true(:),2)^2/numel(x_true);
    end
    
    f1(t+1) = TVnorm(u);
    f2(t+1) = sum(0.5*ysq(:).*exp(-u(:))+u(:));
    obj(t+1) = alpha*TVnorm(u) + sum(0.5*ysq(:).*exp(-u(:))+u(:));
    times(t+1) = cputime - t0;
    %fprintf('DEBUG 4..\n');
%     stopcriterion
    if stopcriterion
        
        switch stopcriterion
            case 1
                criterion = abs( (obj(t+1)-obj(t))/obj(t) );
            case 2
%                 fprintf('debug..')
%                 sum(isnan(clean_prev(:)))
%                 sum(isnan(clean(:)))
%                 norm(clean(:))
%                 norm(clean(:)-clean_prev(:))
%                 criterion = tol + 0.1;
%                 fprintf('you\n')
                criterion = norm(clean(:)-clean_prev(:))/norm(clean_prev(:));
                clean_prev = clean;
            case 3
                criterion = obj(t+1);
            otherwise
                error('Invalid stopping criterion!')
        end
%         criterion
        if criterion < tol
            if verbose
                fprintf('Convergence reached.\n')
            end
            break;
        end
        %fprintf('DEBUG 4a..\n');
    end
    %fprintf('DEBUG 5..\n');
    
    
end
