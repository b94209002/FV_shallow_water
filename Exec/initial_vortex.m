function [u,v] = initial_vortex(mx,my,dx,dy,z)

% subtract the mean
zvec = reshape(z,mx*my,1) - mean(mean(z));

% compute Laplacian matrix in 2D
e =ones(mx,1);
S = spdiags([e e -2*e e e],[-mx+1 -1:1 mx-1],mx,mx);
I = speye(mx);
D = kron(S,I)+kron(I,S); 
D = 1./dx/dy*D;

[L,U]= ilu(D,struct('type','nofill','droptol',1e-02));

pvec = bicgstab(D,zvec,1e-12,500,L,U,zeros(mx*my,1));
psi = reshape(pvec,mx,my);

[v,u] = compute_h_gradient(dx,dy,psi);
v = -v;
return



