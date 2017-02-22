function [sig, iter] = usDenoising_v3(y,varargin)
% -------------------------------------------------------------------------
% USDENOISING: RLTV ULTRASOUND B-MODE FILTERING
% Convex Rayleigh Ultrasound using Log-euclidean Total Variation Priors
% (using embedded C++ code)
% 
% Note: compiled MEX file from C++ source code is used to speed up the
% algorithm
%
% INPUTS: 
% - y:       maximum likelihood estimate
% - alpha:   regularization parameter
% - iMax:    stopping criteria (interrupts program when the algorithm does
%            not converge)
%
% OUTPUTS:
% - sig:     de-noised image
% - speckle: image containing speckle (texture) information
% -------------------------------------------------------------------------

%-----------------   start of: INPUT TESTS ------------------------%
if ~isempty(varargin) && ~isempty(varargin(1)), alpha = varargin{1}; 
else   alpha = 1; end

if length(varargin)>1 && ~isempty(varargin(2)), iMax  = varargin{2}; 
else    iMax = 100; end

if length(varargin)>2 && ~isempty(varargin(3)), tol  = varargin{3}; 
else    tol = 1e-6; end

%-----------------   end of: INPUT TESTS ------------------------%

A = y + eps; %to avoid 0 entries
x = log(A); % monotonic mapping
[l,c,p] = size(x);


% initialize loop
iter = 0;
err = 1;
ev = [];

% First estimate of x:
xs = arrayPrssInit(x, y, alpha);
xs = reshape(xs,l,c,p);

% h = waitbar(0,'Converging to solution...','Name','Denoising');
while (iter < iMax) && (err > tol)
  
    % arrayPrssIterate (estimate_t-1, estimate_t, original_observations, alpha)
    xs = arrayPrssIterate(x, xs, y, alpha);
    xs = reshape(xs,l,c,p);
    diff = x-xs;
    err = norm(diff(:))/(l*c*p);
    x = xs; %assign current estimate to previous one
    %clc;
    iter = iter+1; ev = [ev iter];%disp(ev);
%     waitbar(iter/iMax);

end
% close(h) 

sig = sqrt(exp(x));