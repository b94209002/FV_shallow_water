function sc = BDS_bilinear_interp(c)
% s value on the corner 
% note that s_1/2,1/2 = s(1,1)
%s1 = 49/144;s2 = 7/144;s3=1/144;
s1 = 81/256;s2 = 9/256;s3=1/256;

sc = s1*(c+ circshift(c,[0 1]) +circshift(c,[1 0])+circshift(c,[1 1]));
sc = sc - s2* (circshift(c,[1 2]) +circshift(c,[0 2])+circshift(c,[2 1])+ circshift(c,[-1 1]));
sc = sc - s2* (circshift(c,[2 0]) +circshift(c,[-1 0])+circshift(c,[1 -1])+ circshift(c,[0 -1]));
sc = sc + s3* (circshift(c,[2 2]) +circshift(c,[-1 2])+circshift(c,[2 -1])+ circshift(c,[-1 -1]));
