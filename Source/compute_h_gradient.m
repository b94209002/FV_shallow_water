function [sx sy] = compute_h_gradient(dx,dy,c)

[sxy sx sy sh] = BDS_bilinear_poly(dx,dy,c);

cr = circshift(c,[-1 0]);
cl = circshift(c,[1 0]);

cx = compute_gradient(cr,c,cl)/dx;

sx = -sx;
%sx = -.5*(cr-cl)/dx;

cf = circshift(c,[0 -1]);
cb = circshift(c,[0 1]);

del = .5*(cf-cb)/dy;
%sy = -sign(del).*min(abs(del),abs(sy));
sy = -sy;


