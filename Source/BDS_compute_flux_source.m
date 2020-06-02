function [fxy,fx,fy,fh,div] = BDS_compute_flux_source(u,v,sxy,sx,sy,sh,d)

ix = u>0; iy=v>0;
fxy = BDS_compute_flux_source_loc(ix,iy,sxy);
fx = BDS_compute_flux_source_loc(ix,iy,sx);
fy = BDS_compute_flux_source_loc(ix,iy,sy);
fh = BDS_compute_flux_source_loc(ix,iy,sh);
div = BDS_compute_flux_source_loc(ix,iy,d);






