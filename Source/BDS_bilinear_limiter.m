function [sxy sx sy sh] = BDS_bilinear_limiter(dx,dy,sxy,sx,sy,sh)

% Compute value based on 8a-8d
[LL RL LH RH] = BDS_bilinear_compute_cornor(dx,dy,sxy,sx,sy,sh);

[LL RL LH RH] = BDS_bilinear_constraint(sh,LL,RL,LH,RH);

sx = .5/dx.*((RH+RL)-(LH+LL));
sy = .5/dy.*((LH+RH)-(LL+RL));
sxy = 1./dx./dy.*((RH-RL)-(LH-LL));
return






