function F = BDS_update_2d(dt,dx,dy,umac,vmac,f,c)
%function c = BDS_update_2d(dt,dx,dy,umac,vmac,f,c)


[sxy sx sy sh] = BDS_bilinear_poly(dx,dy,c);

%f0 = zeros(size(f));
[usl usr vsl vsh] = BDS_compute_flux(dt,dx,dy,umac,vmac,sxy,sx,sy,sh,f);

u = umac(1:end-1,:);v = vmac(:,1:end-1);
ua = umac(2:end,:);va = vmac(:,2:end);
%c = c - dt/dx*(ua.*usr-u.*usl) - dt/dy*(va.*vsh - v.*vsl) + f*dt;
F =  - dt/dx.*(ua.*usr-u.*usl) - dt/dy.*(va.*vsh - v.*vsl) ;
return