function U = Heat2(D,f,u0,g1,g2,L,T,nx,nt)
%
% function U = Heat2(D,f,u0,g1,g2,L,T,nx,nt)
%
% The backward difference scheme for solving the initial-boundary value
% problem of the heat equation u_t=D u_{xx}+f for x in [0,L], and
% t in [0,T].
% Input
%      D: the coefficient of the u_{xx} term
%      f=f(x,t): the right hand side function
%      u0=u0(x): the initial value function
%      g1=g1(t): the boundary value function at x=0
%      g2=g2(t): the boundary value function at x=L
%      L:  right end point of the spatial interval
%      T:  right end point of the temporal interval
%      nx: the number of sub-intervals of [0,L]
%      nt: the number of sub-intervals of [0,T]
% Output
%       U: the solution (u_{k,i}), k the time grid point index,
%          i the space grid point index  

% compute some parameters
nx1 = nx+1;
hx = L/nx;
nt1 = nt+1;
ht = T/nt;
r = D*ht/(hx*hx)
r1 = 1+2*r;
%  generate the grid points for x and t variables
xvec = hx*(0:nx);
tvec = ht*(0:nt);
%  generate the entries of the coefficient matrix for the linear 
%  tridiagonal systems to be solved at each time level
avec(2:nx-1) = -r;
bvec(1:nx-1) = r1;
cvec(1:nx-2) = -r;
%  define the dimensions of the unknown 
U = zeros(nt1,nx1);
%  the initial value
U(1,:) = feval(u0,xvec(:)');
%  compute the solution at t=ht; use the program tridiag.m to
%  factor and store A=LU and solve the system at t=ht
U(2,1) = feval(g1,ht);
U(2,nx1) = feval(g2,ht);
fvec(1) = U(1,2)+ht*feval(f,xvec(2),tvec(2))+r*U(2,1);
fvec(2:nx-2) = U(1,3:nx-1)+ht*feval(f,xvec(3:nx-1),tvec(2));
fvec(nx-1) = U(1,nx)+ht*feval(f,xvec(nx),tvec(2))+r*U(2,nx1);
[U(2,2:nx),avec,bvec,ier]=tridiag(avec,bvec,cvec,fvec,nx-1,0);
%  compute the solution at the other time levels
for k = 2:nt1
  U(k,1) = feval(g1,(k-1)*ht);
  U(k,nx1) = feval(g2,(k-1)*ht);
  fvec(1) = U(k-1,2)+ht*feval(f,xvec(2),tvec(k))+r*U(k,1);
  fvec(2:nx-2) = U(k-1,3:nx-1)+ht*feval(f,xvec(3:nx-1),tvec(k));
  fvec(nx-1) = U(k-1,nx)+ht*feval(f,xvec(nx),tvec(k))+r*U(k,nx1);
  U(k,2:nx) = tridiag(avec,bvec,cvec,fvec,nx-1,1);
end
%  plot the numerical solution
surf(xvec,tvec,U)
xlabel('x-axis')
ylabel('t-axis')
zlabel('The numerical solution')
s1 = sprintf('h_t=%6.4f  h_x=%6.4f',ht,hx)
title(['Solution: ',s1])
