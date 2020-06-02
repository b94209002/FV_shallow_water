function c = Strang_update_tracer(mu,mx,my,t,dt,dx,dy,umac,vmac,c,CN,Lr,decay)

% update heating by 1/4 step starting from t
c = update_heating(t,.25*dt,c,Lr,decay);

% diffusion step
% Crank-Nicolson step
%Lc = compute_diffusion(mu,dx,dy,c);
%cs = c + .25*dt*Lc; cs = CN\reshape(cs,mx*my,1);c = reshape(cs,mx,my);

cs = CN\reshape(c,mx*my,1);c = reshape(cs,mx,my);

% update heating by 1/4 step starting from t+.25dt
c = update_heating(t+.25*dt,.25*dt,c,Lr,decay);

% update height by BDS-2d h_t + nabla(uh) = 0
c = c + BDS_update_2d(dt,dx,dy,umac,vmac,zeros(mx,my),c);

% update heating by 1/4 step starting from t+.5dt
c = update_heating(t+.5*dt,.25*dt,c,Lr,decay);

% diffusion step
%Lc = compute_diffusion(mu,dx,dy,c);
%cs = c + .25*dt*Lc; cs = CN\reshape(cs,mx*my,1);c = reshape(cs,mx,my);

cs = CN\reshape(c,mx*my,1);c = reshape(cs,mx,my);

% update heating by 1/4 step starting from t+.75dt
c = update_heating(t+.75*dt,.25*dt,c,Lr,decay);