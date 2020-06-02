function [vuy,uvx] = compute_transverse(v,u,dx,dy,dt) 

% compute vuy
% compute uf and ub
up = circshift(u,[0 -1]);
um = circshift(u,[0 1]);

uy = compute_gradient(up',u',um')'/dy;

%uf_i,j+1/2 and ub_i,j-1/2
uf = u + .5*(dy - dt*v).*uy;uf = [uf(:,end) uf];
ub = u - .5*(dy + dt*v).*uy;ub = [ub ub(:,1)];

%compute vadv_i,j+1/2
vf = circshift(v,[0 -1]);
vb = circshift(v,[0 1]);
vadv = (vb+vf>0 & vb>0).*vb + (vf+vb<0 & vf<0).*vf ; vadv = [vadv(:,end) vadv];

uhat = (vadv > 0).*ub + (vadv < 0).*uf + .5*(vadv == 0).*(uf + ub);

odx = .5/dx;
vuy = odx*(vadv(:,2:end)+vadv(:,1:end-1)).*(uhat(:,2:end) - uhat(:,1:end-1));


% compute uvx
% compute vl and vr
vp = circshift(v,[-1 0]);
vm = circshift(v,[1 0]);

vx = compute_gradient(vp,v,vm)/dx;

vl = v + .5*(dx - dt*u).*vx;vl = [vl(end,:); vl];
vr = v - .5*(dx + dt*u).*vx;vr = [vr; vr(1,:)];

% compute uadv 
ul = circshift(u,[-1 0]);
ur = circshift(u,[1 0]);
uadv = (ur+ul>0 & ur>0).*ur + (ur+ul<0 & ul<0).*ul ; uadv = [uadv(end,:); uadv];

vhat = (uadv > 0).*vr + (uadv < 0).*vl + .5*(uadv == 0).*(vr + vl);

uvx = odx*(uadv(2:end,:)+uadv(1:end-1,:)).*(vhat(2:end,:) - vhat(1:end-1,:));
