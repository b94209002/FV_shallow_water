function [LL RL LH RH] = BDS_bilinear_compute_cornor(dx,dy,sxy,sx,sy,sh)
% compute corner value based on bilinear polynomial  
LL = sh - .5*dx.*sx -.5*dy.*sy + .25*dx.*dy.*sxy; 
LH = sh - .5*dx.*sx +.5*dy.*sy - .25*dx.*dy.*sxy; 
RL = sh + .5*dx.*sx -.5*dy.*sy - .25*dx.*dy.*sxy; 
RH = sh + .5*dx.*sx +.5*dy.*sy + .25*dx.*dy.*sxy; 