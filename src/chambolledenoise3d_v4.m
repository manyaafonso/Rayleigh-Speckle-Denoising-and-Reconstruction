function fe = chambolledenoise3d_v4(g, alpha, iter )

[M,N,L] = size(g);

p = zeros(M,N,L,3);

tau = 0.1;

for k = 1:iter
    
    %div_p = divergence3D(p(:,:,:,1),p(:,:,:,2),p(:,:,:,3));
    div_p = divergence3D_v2(p(:,:,:,1),p(:,:,:,2),p(:,:,:,3));
    
    [t1,t2,t3] = gradient3D_v2(div_p - (1/alpha)*g);
    tmp_grad(:,:,:,1) = t1;
    tmp_grad(:,:,:,2) = t2;
    tmp_grad(:,:,:,3) = t3;
    
    num = p + tau*tmp_grad;
    
    den = 1 + tau*sqrt(tmp_grad(:,:,:,1).^2 + tmp_grad(:,:,:,2).^2 + tmp_grad(:,:,:,3).^2 );
    
    for i = 1:3
        p(:,:,:,i) = num(:,:,:,i)./den;
    end
    
end

fe = g - alpha*divergence3D_v2(p(:,:,:,1),p(:,:,:,2),p(:,:,:,3));
