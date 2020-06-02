function c = update_heating(t,dt,c,Lr,decay)
% apply Runge-Kutta 2nd for nonautonomous equations to update heating
% page 125 of Randall J. Leveque Finite Difference method
% Lr       : L in Hendricks et al. 2014
% decay(t) : Q0*alpha*t*exp(alpha*t) in Hendricks et al. 2014

% first step of RK
Q = Lr * decay(t);
ch = c + .5*dt*Q.*c;

% final step of RK
Q = Lr * decay(t+.5*dt);
c = c + dt*Q.*ch;