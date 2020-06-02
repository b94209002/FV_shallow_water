function [u,v] = CN_update_momentum(m,dt,dx,dy,u,v,umac,vmac,Lu,Lv,Sx,Sy,CN)

us = u + BDS_update_2d(dt,dx,dy,umac,vmac,Sx+Lu,u) + dt*( Sx + .5*Lu);
vs = v + BDS_update_2d(dt,dx,dy,umac,vmac,Sy+Lv,v) + dt*( Sy + .5*Lv);

us = CN\reshape(us,m*m,1);u = reshape(us,m,m);
vs = CN\reshape(vs,m*m,1);v = reshape(vs,m,m);