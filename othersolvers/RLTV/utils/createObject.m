function v = createObject(n,bg,fg,mode,noise)
% v = createObject(n,bg,fg,mode) creates simulated ultrasound volume with
% a specific object given in 'mode'.
%
%  Input:   n - Matrix [1 x3] with element size
%           bg - background intensity;
%           fg - foreground intensity;
%           mode - 'cube', 'tube', 'sphere'
%           noise - binary entry if Rayleigh noise is added or not
%           
%  Output: 3D matrix array with slices intensities

%-----------------   start of: INPUT TESTS ------------------------%
if ~exist('n','var') || isempty(n), n = 64; end
if ~exist('bg','var') || isempty(bg), bg = 10; end
if ~exist('fg','var') || isempty(fg), fg = 100; end
if ~exist('mode','var') || isempty(mode), mode = 'tube'; end
if ~exist('noise','var') || isempty(noise), noise = 1; end
%-----------------   end of: INPUT TESTS ------------------------%

v = ones(n(1),n(2),n(3))*bg;
switch mode,
    case 'cube'
        v(n(1)/2-n(1)/4:n(1)/2+n(1)/4, ...
            n(2)/2-n(2)/4:(2)/2+n(2)/4, ...
            n(3)/2-n(3)/4:n(3)/2+n(3)/4) = fg;
             
    case 'tube'
        [x,y,z] = meshgrid(linspace(-1,1,n(1)),linspace(-1,1,n(2)),linspace(-1,1,n(3)));
        v = double(sqrt(x.^2+y.^2) <= 0.75 & sqrt(x.^2+y.^2) >= 0.35).*fg +...
            double(sqrt(x.^2+y.^2) > 0.75 & sqrt(x.^2+y.^2) < 0.35).*bg;
                
    case 'sphere'
        [x,y,z] = meshgrid(linspace(-1,1,n(1)),linspace(-1,1,n(2)),linspace(-1,1,n(3)));
        v = double(sqrt(x.^2+y.^2+z.^2) <= 0.75).*fg+...
            double(sqrt(x.^2+y.^2+z.^2) > 0.75).*bg;
        
    otherwise
        errordlg('Not a valid object');
        return;
end
if noise, v = raylrnd(v); end %Random Rayleigh noise corruption
