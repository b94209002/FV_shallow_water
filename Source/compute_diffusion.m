function D = compute_diffusion(mu,dx,dy,c)

D = circshift(c,[-1,0])+circshift(c,[1,0]);
D = D+circshift(c,[0,-1])+circshift(c,[0,1]);
D = D-4*c;
D = mu*D/dx/dy;