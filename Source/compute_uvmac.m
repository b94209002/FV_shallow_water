function [umac,vmac] = compute_uvmac(dt,dx,dy,u,v,Sx,Sy)

[vuy,uvx]= compute_transverse(v,u,dx,dy,dt);

up = circshift(u,[-1 0]);
um = circshift(u,[1 0]);

ux = compute_gradient(up,u,um)/dx;

ul = u - .5*ux.*(dx - u*dt) -.5*dt*vuy + .5*dt*Sx;ul = [ul(end,:); ul];
ur = u - .5*ux.*(dx + u*dt) -.5*dt*vuy + .5*dt*Sx;ur = [ur; ur(1,:)];

umac = (ul+ur>0 & ul>0).*ul + (ul+ur<0 & ur<0).*ur ;

vp = circshift(v,[0 -1]);
vm = circshift(v,[0 1]);

vy = compute_gradient(vp',v',vm')'/dy;

vt = v - .5*vy.*(dy - v*dt) -.5*dt*uvx + .5*dt*Sy;vt = [vt(:,end) vt];
vb = v - .5*vy.*(dy + v*dt) -.5*dt*uvx + .5*dt*Sy;vb = [vb vb(:,1)];

vmac = (vt+vb>0 & vt>0).*vt + (vt+vb<0 & vb<0).*vb ;