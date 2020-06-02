function [temp,kdp,sumdiff] = BDS_bilinear_iter(temp,kdp,sgn,sumdiff,mn,mx)
% determine the max/min by the sign of sumdiff
%extre = .5*((sign(sumdiff) +1)*mn-(sign(sumdiff)-1)*mx); 
%redfac = min( sumdiff/kdp, temp - extre);
tmp =.5*(sgn+1)*(temp - mn) - .5*(sgn-1)*(mx - temp); 
redfac = sgn*sumdiff/min(1,kdp);
redfac = min(redfac,tmp);
kdp = kdp - 1;
sumdiff = sumdiff - sgn*redfac;
temp = temp - sgn*redfac;
    
