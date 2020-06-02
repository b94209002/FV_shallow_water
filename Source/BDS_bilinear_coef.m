function [sxy sx sy] = BDS_bilinear_coef(dx,dy,sc)

sx = .5/dx.*(circshift(sc,[-1,-1])+circshift(sc,[-1,0])-(circshift(sc,[0,-1])+circshift(sc,[0 0])));
sy = .5/dy.*(circshift(sc,[-1,-1])-circshift(sc,[-1,0])+(circshift(sc,[0,-1])-circshift(sc,[0,0])));
sxy = 1/dx./dy.*((circshift(sc,[-1,-1])-circshift(sc,[-1,0]))-circshift(sc,[0,-1])+circshift(sc,[0,0]));