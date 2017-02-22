function [clean, obj, times, mses] = rayleighReconstruct2D_v2(varargin)


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
    maxiters = 8;
    chambolleit = 5;
    inner_iters = 4;
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

[M,N] = size(observed);

maskvec = mask(:);

A = @(x)  x(maskvec==1);
AT = @(x) padwithzeros(x,maskvec);%x(:,:,q);

vox = observed(:);
y = vox(maskvec==1);

ysq = y.^2;

% x = zeros(L,M,N);
% x = x(:);
x = 0.5*vox.^2;

Ax = A(x);

u = Ax;
expu = exp(u);
du = u;
v = 0*x;
dv = v;

H(1) = sum(sum(0.5*ysq.*exp(-u)+u));
obj(1) = alpha*TVnorm(reshape(v,M,N)) + H(1);
    
times(1) = 0;
t0 = cputime;
mses = [];

clean_prev = x;

kappa = mu2/mu1;
invLS = (1/(1+kappa))*maskvec + (1/kappa)*(1-maskvec);

for t = 1:maxiters
    if verbose
        t
    end
    
    r = Ax-du;
    
%     for i = 1:length(y)
%         u(i) = 0;
%         for nm = 1:inner_iters
%             u(i) = u(i) - (-0.5*ysq(i)*exp(-u(i))+1-mu1*exp(u(i))*(r(i)-exp(u(i))))/( 0.5*ysq(i)*exp(-u(i)) - mu1*r(i)*exp(u(i))+2*mu1*exp(2*u(i)) );
%         end
%     end

    u = 0*u;
    expu = exp(u);
    for nm = 1:inner_iters
        %u = u - ( -0.5*ysq.*exp(-u)+1-mu1*exp(u).*(r-exp(u)) )./( 0.5*ysq.*exp(-u) - mu1*r.*exp(u)+2*mu1*exp(2*u) );
        u = u - ( -0.5*ysq./expu+1-mu1*expu.*(r-expu) )./( 0.5*ysq./expu - mu1*r.*expu+2*mu1*(expu.^2) );
        expu = exp(u);
    end
    
    v = projk(reshape(x-dv,M,N), alpha/mu2, chambolleit );
    
%     alpha = alpha/1.02;
    
    v = v(:);
    
    %s = AT(exp(u)+du) + (mu2/mu1)*(v+dv);
    s = AT(expu+du) + (mu2/mu1)*(v+dv);
    %x = s -(1/(1+mu2/mu1))*AT(A(s));
    x = invLS.*s;
%     x = max(x,0);
    
    Ax = A(x);
    
    %du = du - (Ax-exp(u));
    du = du - (Ax-expu);
    dv = dv - (x-v);
    
    clean = reshape(x,M,N);
    
    if compute_mse
        mses(t) = norm(clean(:)-x_true(:),2)^2/numel(x_true);
    end
    
    %H(t+1) = sum(sum(sum(0.5*ysq.*exp(-u)+u)));
    H(t+1) = sum(sum(0.5*ysq./(Ax+1e-6)+log(Ax+1e-6)));
    %obj(t+1) = alpha*TVnorm3d(v) + H(t+1);
%     obj(
%     if ~isreal(obj(t+1))
%         isreal(H(t+1))
%         sum(isreal(x(:)))
%         sum(isreal(u(:)))
%         sum(isreal(v(:)))
%         isreal(0.5*ysq./(A(x)+1e-6))
%         isreal( log(A(x)+1e-3) )
%         break
%     end
%     Ax = A(x);
%     eu = exp(u);
    %
    obj(t+1) = alpha*TVnorm(reshape(x,M,N)) + H(t+1);% + mu1*norm(Ax(:)-eu(:)-du(:))^2 + mu2*norm(x(:)-v(:)-dv(:))^2;
    
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
