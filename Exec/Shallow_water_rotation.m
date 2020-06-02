clear
addpath('../source')
exact_sol =1;
%exact_sol = 0;mx = fliplr([100 200 400 800 1600]/2); hx= 1./mx;

ax = -5; bx = 5; ay = ax; by = bx;
m = 100;mx = m; my = m; % number of grid points
dx = (bx-ax)/mx; dy = (by - ay)/my; % Grid size
x = linspace(ax+dx/2,bx-dx/2,m)';y = linspace(ay+dy/2,by-dy/2,m)';
[xx yy] = meshgrid(x,y); %center

%xi = linspace(h/2,1-h/2,m)';% x-coordinate on the interface
g = 0.1; % advection coiefficient
H = 10;
L = bx-ax;
nu = 0.05; 
mu = 0.1;
f = 0.01;
dt = nu*dx;
T = 1;
alpha = 1; % cooling rate
Q0 = -0.5; %initial heating
xc = 0.; yc = 0.; % center of the vortex;
r1 = 0.45; r2 = 0.55; r3 = 0.95; r4 = 1.05;
z1 = 1e-1; z2 = 0;

% initial condition for u,v and c 
%[u,v] = fixed_velocity(mx,my,0,0);
r = sqrt((xx-xc).^2 + (yy-yc).^2);
S = @(x) 1 - 3*x.^2 + 2*x.^3;
zeta = z1*(r<=r3) + z2*(r>=r4);
zeta = zeta + (r>r3 & r<r4).*(z1*S((r-r3)/(r4-r3))+ z2*S((r4-r)/r4-r3));
[u,v] = initial_vortex(mx,my,dx,dy,zeta);
c = H*ones(mx,my);% + .5*ones(mx,my).*(sqrt((xx-.5).^2+(yy-.5).^2)<=0.1);

% computing heating 
decay = @(t) Q0*alpha*t.*exp(-alpha*t);
Lr = (r <= r1 | r >= r4)*0 + (r >= r2 & r <= r3)*1 ;
Lr = Lr + (r1 < r & r < r2).*S((r2-r)/(r2-r1)) + (r3 < r & r < r4).*S((r-r3)/(r4-r3));

CN = Crank_Nicolson_A(m,mu,dt,dx,dy,-.5);
BE = Crank_Nicolson_A(m,mu,.5*dt,dx,dy,-1);
% compute the height gradient at initial step
[sx,sy] = compute_h_gradient(dx,dy,c);
i=0;uh = zeros(m,1);hh = uh; vh = zeros(m+1,1);sh = uh;


for t =dt:dt:T

Lu = compute_diffusion(mu,dx,dy,u);
Lv = compute_diffusion(mu,dx,dy,v);
Fx = g*sx+Lu+f*v; 
Fy = g*sy+Lv-f*u;

% compute umac, vmac
[umac,vmac] = compute_uvmac(dt,dx,dy,u,v,Fx,Fy);
% compute velocity at grids and t = n+1/2
[um,vm] = compute_velocity_middle(umac,vmac);
Div = (umac(2:end,:) - umac(1:end-1,:))/dx + (vmac(:,2:end) - vmac(:,1:end-1))/dy;


% update height by Strang splitting with diffusion and advection.
c = Strang_update_tracer(mu,mx,my,t,dt,dx,dy,umac,vmac,c,BE,Lr,decay);

% compute height gradient again at t = n + 1
osx = sx; osy = sy;
[sx,sy] = compute_h_gradient(dx,dy,c);

Fx = .5*g*(osx+sx) + f*vm;
Fy = .5*g*(osy+sy) - f*um;

[u,v] = CN_update_momentum(m,dt,dx,dy,u,v,umac,vmac,Lu,Lv,Fx,Fy,CN);

subplot(221);pcolor(xx,yy,c'-H);shading flat;colorbar;
title([' t = ' num2str(t) ' max = ' num2str(max(max(c-H))) ', min = ' num2str(min(min(c-H))) ] );
subplot(222);pcolor(xx,yy,u');shading flat;colorbar;title('u');%caxis([-1 1]);
subplot(223);pcolor(xx,yy,v');shading flat;colorbar;title('v');%caxis([-1 1]);
subplot(224);pcolor(xx,yy,Div');shading flat;colorbar;title('\delta')
drawnow
%pcolor(xx,yy,c'-H);title(num2str(t));caxis([-1 1]);colorbar;shading flat
%print(['fig/time_step' num2str(t*1000,'%04i')],'-dpng','-r300')
end

%compute_norm(wnorm1,wnorm2,wnorm3,wnorm4,wnorm5)

return
