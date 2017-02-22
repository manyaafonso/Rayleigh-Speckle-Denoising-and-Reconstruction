[x,y,z] = meshgrid(linspace(-1,1,n),linspace(-1,1,n),linspace(-1,1,n));
figure('name','Image selection');
subplot(121);
p=patch(isosurface(x,y,z,v,1));
set(p,'facecolor','red','edgecolor','none');
alpha(0.5);
box on
hold on
axis on
rotate3d on
aux(cg<15)=NaN;
i=19;  mesh(cx(:,:,i),cy(:,:,i),cz(:,:,i),cg(:,:,i)), grid off; axis tight; axis square
view(-5,10)

subplot(122); imshow(cg(:,:,i),[]);