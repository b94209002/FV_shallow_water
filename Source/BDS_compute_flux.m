function [usl usr vsl vsh] = BDS_compute_flux(dt,dx,dy,umac,vmac,sxy,sx,sy,sh,f)
u = umac(1:end-1,:);v = vmac(:,1:end-1);
ua = umac(2:end,:);va = vmac(:,2:end);
% ua va are velocity at i+1/2 and j+1/2
%ua = circshift(u,[-1 0]);va = circshift(v,[0 -1]);
div = (ua-u)./dx+ (va-v)./dy;
% corner index dx/2 or -dx/2 us, vs
us = 2*(ua>0)-1;
vs = 2*(va>0)-1;

%usx = ((ua>0).*sl + (ua<=0).*sr);
% DEF area in BDS paper
% u at u(i+1/2) vs  v(j+1/2)
vv = .5*((us+1).*va + (1-us).*circshift(va,[-1,0]));
[fxy,fx,fy,fh,fdiv] = BDS_compute_flux_source(ua,vv,sxy,sx,sy,sh,div);
%[fxy,fx,fy,fh,fdiv] = BDS_compute_flux_source(ua,vv,sxy,sx,sy,sh,f);
vt = 2*(vv>0) -1;
up = .5*((1+vt).*ua + .5*(1-vt).*(1+us.*sign(circshift(us,[0,-1]))).*circshift(ua,[0,-1]));
Gammayp = fh + 1/6*(fx.*(3*us*dx-2*dt*(ua+up)) + fy.*(3*vt*dy-2*dt*vv));
Gammayp = Gammayp + 1/12*fxy.*(3*dx*dy*us.*vt - 2*dt*dy*vt.*(ua+up) - 2*dt*dx*us.*vv + dt*dt*(ua+2*up).*vv);
Gammayp = Gammayp .* (1 - dt/3*fdiv); 
%Gammayp = Gammayp + dt/3 * fdiv;

% u at u(i+1/2) vs  v(j-1/2)
vv = .5*((us+1).*v + (1-us).*circshift(v,[-1,0]));
[fxy,fx,fy,fh,fdiv] = BDS_compute_flux_source(ua,vv,circshift(sxy,[0 1]),circshift(sx,[0 1]),circshift(sy,[0 1]),circshift(sh,[0 1]),circshift(div,[0 1]));
%[fxy,fx,fy,fh,fdiv] = BDS_compute_flux_source(ua,vv,circshift(sxy,[0 1]),circshift(sx,[0 1]),circshift(sy,[0 1]),circshift(sh,[0 1]),circshift(f,[0 1]));
vt = 2*(vv>0)-1;
up = .5*((1-vt).*ua + .5*(1+vt).*(1+us.*sign(circshift(us,[0,1]))).*circshift(ua,[0,1]));
Gammaym = fh + 1/6*(fx.*(3*us*dx-2*dt*(ua+up)) + fy.*(3*vt*dy-2*dt*vv));
Gammaym = Gammaym + 1/12*fxy.*(3*dx*dy*us.*vt - 2*dt*dy*vt.*(ua+up) - 2*dt*dx*us.*vv + dt*dt*(ua+2*up).*vv);
Gammaym = Gammaym .* (1 - dt/3*fdiv); 
%Gammaym = Gammaym + dt/3*fdiv; 

%note here suppose u> 0
sxp = .5*sx.*(dx-ua*dt) + sh;
sxm = .5*circshift(sx,[-1,0]).*(-dx-ua*dt) + circshift(sh,[-1,0]);
%First two term of 17
sl = sxp .* ( 1-.5*dt/dx.*(ua - u)) + .5*dt*f; 
sr = sxm .* ( 1-.5*dt/dx*(circshift(ua,[-1,0]) - ua)) +.5*dt*circshift(f,[-1,0]);

usx = .5*((us+1).*sl + (1-us).*sr);

vv = .5*((us+1).*va + (1-us).*circshift(va,[-1,0]));
vm = .5*((us+1).*v + (1-us).*circshift(v,[-1,0]));
%circshift(vv,[0 1]);
uslr = .5*dt/dy*(vv.*Gammayp - vm.*Gammaym);

usr = ( usx - uslr );
usl = circshift(usr,[1,0]);

%note here suppose v > 0 
syp = .5*sy.*(dy-va*dt) + sh;
sym = .5*circshift(sy,[0,-1]).*(-dy-va*dt) + circshift(sh,[0,-1]);
st = syp .* ( 1-.5*dt/dy.*(va - v)) +.5*dt*f;
sb = sym .* ( 1-.5*dt/dy*(circshift(va,[0 -1]) - va)) + .5*dt*circshift(f,[0 -1]);

vsx = .5*((vs+1).*st + (1- vs).*sb);

uu = .5*((vs+1).*ua + (1-vs).*circshift(ua,[0,-1]));
[fxy,fx,fy,fh,fdiv] = BDS_compute_flux_source(uu,va,sxy,sx,sy,sh,div);
%[fxy,fx,fy,fh,fdiv] = BDS_compute_flux_source(uu,va,sxy,sx,sy,sh,f);
ut = 2*(uu>0) -1;
vp = .5*((1+ut).*va + .5*(1-ut).*(1+vs.*sign(circshift(vs,[-1,0]))).*circshift(va,[-1,0]));
Gammaxp = fh + 1/6*(fy.*(3*vs*dy-2*dt*(va+vp)) + fx.*(3*ut*dx-2*dt*uu));
Gammaxp = Gammaxp + 1/12*fxy.*(3*ut.*vs*dx*dy - 2*dt*dx*ut.*(va+vp) - 2*dt*dy*vs.*uu + dt*dt*uu.*(va+2*vp));
Gammaxp = Gammaxp .* (1-dt/3*fdiv);
%Gammaxp = Gammaxp + dt/3*fdiv;

uu = .5*((vs+1).*u + (1-vs).*circshift(u,[0,-1]));
[fxy,fx,fy,fh,fdiv] = BDS_compute_flux_source(uu,va,circshift(sxy,[1 0]),circshift(sx,[1 0]),circshift(sy,[1 0]),circshift(sh,[1 0]),circshift(div,[1 0]));
%[fxy,fx,fy,fh,fdiv] = BDS_compute_flux_source(uu,va,circshift(sxy,[1 0]),circshift(sx,[1 0]),circshift(sy,[1 0]),circshift(sh,[1 0]),circshift(f,[1 0]));
ut = 2*(uu>0) -1;
vp = .5*((1-ut).*va + .5*(1+ut).*(1+vs.*sign(circshift(vs,[1,0]))).*circshift(va,[1,0]));
Gammaxm = fh + 1/6*(fy.*(3*vs*dy-2*dt*(va+vp)) + fx.*(3*ut*dx-2*dt*uu));
Gammaxm = Gammaxm + 1/12*fxy.*(3*ut.*vs*dx*dy - 2*dt*dx*ut.*(va+vp) - 2*dt*dy*vs.*uu + dt*dt*uu.*(va+2*vp));
Gammaxm = Gammaxm .* (1-dt/3*fdiv);
%Gammaxm = Gammaxm + dt/3*fdiv;

uu = .5*((vs+1).*ua + (1-vs).*circshift(ua,[0,-1]));
um = .5*((vs+1).*u + (1-vs).*circshift(u,[0,-1]));
vshx = .5*dt/dx.*(uu.*Gammaxp - um.*Gammaxm);

vsh = ( vsx - vshx);
vsl = circshift(vsh,[0,1]);
          
return

