function planes = arbitrarySlices(noPlanes,v,n,op)
% generates ultrasound slices from a volume of interest by sweeping along
% it with randomly rotation
%
% inputs:   noPlanes - number of random planes to acquire
%           v        - original volume of data
%           n        - number of grid elements per dimension
%           op       - 'random', 'regular' or 'jitter' image acquisition
%
% outputs:   planes - structure containing the coordinates and grey levels of each
% acquired plane

%-----------------   start of: INPUT TESTS ------------------------%
if ~exist('n','var') || isempty(n), n = 64; end
if ~exist('noPlanes','var') || isempty(noPlanes), noPlanes = 1; end
if ~exist('op','var') || isempty(op), op = 'jitter'; end
%-----------------   end of: INPUT TESTS ------------------------%

switch op
    case 'random',
        rval = -90:5:90;
        rx = rval(randperm(numel(rval)));
        ry = rval(randperm(numel(rval)));
        rz = rval(randperm(numel(rval)));
        
    case 'regular',
        rx = 10*ones(noPlanes,1);
        ry = 30*ones(noPlanes,1);
        rz = 10*ones(noPlanes,1);
        
    case 'jitter'
        amount = 2;
        rx =  10*ones(noPlanes,1)+ amount*randn(noPlanes,1);
        ry = -10*ones(noPlanes,1)+ amount*randn(noPlanes,1);
        rz =  10*ones(noPlanes,1)+ amount*randn(noPlanes,1);
        
    otherwise,
        errordlg('Not a valid mode');
        return;
end


planes = []; % structure which will contain the coordinates and grey values of each acquired plane


center = [0 0 0]; % center of rotation
xdir = [1 0 0]; ydir = [0 1 0]; zdir = [0 0 1]; % rotate along these directions

% volume coordinates
[x,y,z] = meshgrid(linspace(-1,1,n),linspace(-1,1,n),linspace(-1,1,n));

figure('OuterPosition',get(0, 'ScreenSize'));

idx = 0;
for i = linspace(-1,1,noPlanes),

    idx = idx + 1;
    % creates a default surface and rotates it 
    hsp = surf(linspace(-1,1,n),linspace(-1,1,n),zeros(n)+i);    %constant speed Transladation%
    
    rotate(hsp,xdir,rx(idx),center)
    rotate(hsp,ydir,ry(idx),center)
    rotate(hsp,zdir,rz(idx),center)

    xd = get(hsp,'XData'); % plane coordinates
    yd = get(hsp,'YData');
    zd = get(hsp,'ZData');
    delete(hsp)
    
    subplot(121), 
    h = slice(x,y,z,v,xd,yd,zd);
    title('observed volume');
    pd = get(h,'CData'); % plane grey values
    
    %  ----  Object surface representation  -----%
    p=patch(isosurface(x,y,z,v,1));
    set(p,'facecolor','red','edgecolor','none');
    alpha(0.5);
%     box on;
%     grid on;

    %------------------------%

    hold off
    axis([-1.5 1.5 -1.5 1.5 -1.5 1.5]); axis square
    colormap gray
    view(-5,10)
%     axis off;
    
    subplot(122), 
    imshow(pd,[]), title('slice');
    drawnow;
    
    planes(idx).x = xd;
    planes(idx).y = yd;
    planes(idx).z = zd;
    planes(idx).g = pd;    
    
end