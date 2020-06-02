function D = Crank_Nicolson_A(m,mu,dt,dx,dy,c)

%note the vector in spdiags will be shifted.
e =ones(m,1);
S = spdiags([e e -2*e e e],[-m+1 -1:1 m-1],m,m);
I = speye(m);
D = kron(S,I)+kron(I,S); 
D = c*mu*dt/dx/dy*D;% c=.5 for forward; c=-.5 for backward
D = speye(m*m) + D;