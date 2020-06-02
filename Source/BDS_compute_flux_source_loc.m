function f = BDS_compute_flux_source_loc(ix,iy,s)

f = ix.*iy.*s;
f = f + (1-ix).*iy.*circshift(s,[-1,0]);
f = f + ix.*(1-iy).*circshift(s,[0,-1]);
f = f + (1-ix).*(1-iy).*circshift(s,[-1,-1]);