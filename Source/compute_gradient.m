function ux = compute_gradient(up,u,um)

del = .5*(up-um);
dpls = 2*(up-u);
dmin = 2*(u-um);

slim = min(abs(dpls),abs(dmin));
slim = (dpls.*dmin > 0).*slim;

ux = sign(del).*min(slim,abs(del));
