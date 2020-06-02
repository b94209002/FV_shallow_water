function [sxy sx sy sh] = BDS_bilinear_poly(dx,dy,c)

sc = BDS_bilinear_interp(c);

[sxy sx sy] = BDS_bilinear_coef(dx,dy,sc);
sh = c;

[sxy sx sy sh] = BDS_bilinear_limiter(dx,dy,sxy,sx,sy,sh);
return


