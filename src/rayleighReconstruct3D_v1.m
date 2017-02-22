function [clean, obj, times, mses] = rayleighReconstruct3D_v1(varargin)


if nargin < 5
    error('Must specify observed volume y, mask, regularization parameter alpha, and weights mu1, mu2!')
else
    observed = varargin{1};
    mask = varargin{2};
    alpha = varargin{3};
    mu1 = varargin{4};
    mu2 = varargin{5};
    
    % default parameters
    compute_mse = 0;
    maxiters = 200;
    chambolleit = 5;
    inner_iters = 5;
    verbose = 0;
    stopcriterion = 0;

    if nargin == 6
        %%% if options are provided, use them, else use default parameters.
        opts = varargin{6};
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

[L,M,N] = size(observed);

maskvec = mask(:);

A = @(x)  x(maskvec==1);
AT = @(x) padwithzeros(x,maskvec);%x(:,:,q);

vox = observed(:);
y = vox(maskvec==1);

ysq = y.^2;

x = zeros(L,M,N);
x = x(:);

u = A(x);
du = u;
v = x;
dv = v;

H(1) = sum(sum(sum(0.5*ysq.*exp(-u)+u)));
obj(1) = alpha*TVnorm3d(x) + H(1);
    
times(1) = 0;
t0 = cputime;
mses = [];

clean_prev = x;

for t = 1:maxiters
    if verbose
        t
    end
    
    r = A(x)-du;
    for i = 1:length(y)
        u(i) = 0;
        for nm = 1:inner_iters
            u(i) = u(i) - (-0.5*ysq(i)*exp(-u(i))+1-mu1*exp(u(i))*(r(i)-exp(u(i))))/( 0.5*ysq(i)*exp(-u(i)) - mu1*r(i)*exp(u(i))+2*mu1*exp(2*u(i)) );
        end
    end
    
    %v = chambolledenoise3d_v2(reshape(x-dv,M,N,L), alpha/mu2, chambolleit );
    v = chambolledenoise3d_v4(reshape(x-dv,M,N,L), alpha/mu2, chambolleit );
    v = v(:);
    
    s = AT(exp(u)+du) + (mu2/mu1)*(v+dv);
    x = s -(1/(1+mu2/mu1))*AT(A(s));
    
    du = du - (A(x)-exp(u));
    dv = dv - (x-v);
    
    clean = x;
    
    if compute_mse
        mses(t) = norm(clean(:)-x_true(:),2)^2/numel(x_true);
    end
    
    H(t+1) = sum(sum(sum(0.5*ysq.*exp(-u)+u)));
    obj(t+1) = alpha*TVnorm3d(x) + H(t);
    
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
