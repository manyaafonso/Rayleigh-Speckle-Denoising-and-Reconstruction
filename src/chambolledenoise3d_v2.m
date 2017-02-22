function fe = chambolledenoise3d_v2(g, alpha, iter )

[M,N,L] = size(g);

p = zeros(M,N,L,3);

tau = 0.1;

for k = 1:iter
    
%     tic
    div_p = divergence3d(p(:,:,:,1),p(:,:,:,2),p(:,:,:,3));
%     toc
%     tic
    tmp_grad = gradient_3d(div_p - (1/alpha)*g);
%     toc
    
    num = p + tau*tmp_grad;
    
%     tic

    den = 1 + tau*sqrt(tmp_grad(:,:,:,1).^2 + tmp_grad(:,:,:,2).^2 + tmp_grad(:,:,:,3).^2 );
%     toc
%     
%     tic

    for i = 1:3
        p(:,:,:,i) = num(:,:,:,i)./den;
    end
%     toc
    
end

fe = g - alpha*divergence3d(p(:,:,:,1),p(:,:,:,2),p(:,:,:,3));
