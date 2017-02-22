function view3Dmodified(seg)

visSeg = seg;
visSeg = permute(visSeg, [2 1 3]);
%visSeg = flipdim(visSeg, 1);
%visSeg = flipdim(visSeg, 2);
flippedx = inline(sprintf('%d - x + 1', size(visSeg, 2))); 
%The 2nd dimension of the matrix ends up being 'called' x by matlab.
%wtf?

sizeSeg = size(visSeg);
sn = round(sizeSeg/2);
% sn(1)
% sn(2)
% sn(3)
hslc=slice(visSeg, flippedx(sn(1)) ,sn(2),sn(3));
axis equal; axis vis3d; axis off;
colormap gray;
set(hslc(1:3),'LineStyle','none');
%xlabel 'p-a' ;ylabel 'l-r' ;zlabel 'i-s';
xlabel 'x' ;ylabel 'y' ;zlabel 'z';

%set(gca, 'XDir', 'reverse'); 
%The above XDir didn't work, without screwing up the data.  So, I'll
%just change the labels...WAIT!: that's wrong because the slice numbers
%become wrong then.  Have to flip...
myx = get(gca, 'XTickLabel');
set(gca, 'XTickLabel', flippedx(str2num(myx)));

%title 'Manually enable/disable Rotate3D.'
%To maintain the 3D perspective throughout.
camPos3D = gcf;%get(gcf, 'CameraPosition');


%Too many problems trying to regain control from the rotate3d.
%rotHandle = rotate3d(handles{4}); so it's up to the user.
rotHandle = rotate3d;
% setAllowAxesRotate(rotHandle, handles{1}, false);
% setAllowAxesRotate(rotHandle, handles{2}, false);
% setAllowAxesRotate(rotHandle, handles{3}, false);
%End of the 3D view elements
