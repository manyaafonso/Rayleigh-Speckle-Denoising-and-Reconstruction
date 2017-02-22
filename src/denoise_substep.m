function u = denoise_substep(y, mu, r, inner_iters)

u = zeros(size(y));%log(y+1e-3);
ysq = y.^2;

for i = 1:inner_iters
    expmu = exp(-u);
    u = u - (mu*(u-r)+1-0.5*ysq.*expmu)./(mu+0.5*ysq.*expmu);
end
