function U = Wave(a,f,u0,v0,g1,g2,L,T,nx,nt)
%
% function U=Wave(a,f,u0,v0,g1,g2,L,T,nx,nt)
%
% The centered difference scheme for solving the initial-boundary value
% problem of the wave equation u_{tt}=a u_{xx}+f for x in [0,L], and
% t in [0,T].
% Input
%      a: the coefficient of the u_{xx} term
%      f=f(x,t): the right hand side function
%      u0=u0(x): the initial value function
%      v0=v0(x): the initial derivative value function
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
ht2 = ht*ht;
r = a*ht*ht/(hx*hx)
r1 = 2*(1-r);
if r > 1
  disp('The ratio is too big; the scheme is unstable!')
end
%  define the grid points
xvec = hx*(0:nx);
tvec = ht*(0:nt);
%  compute the solution at t=0 and ht using the initial values
U = zeros(nt1,nx1);
U(1,:) = feval(u0,xvec);
U(2,1) = feval(g1,ht);
U(2,2:nx) = 0.5*(r*(U(1,1:nx-1)+U(1,3:nx+1))+r1*U(1,2:nx) ...
            +ht2*feval(f,xvec(2:nx),tvec(1))) ...
            +ht*feval(v0,xvec(2:nx));
U(2,nx1) = feval(g2,ht);
%  compute the solution at the other times
for k = 2:nt
  U(k+1,1) = feval(g1,k*ht);
  U(k+1,nx1) = feval(g2,k*ht);
  U(k+1,2:nx) = r*(U(k,1:nx-1)+U(k,3:nx+1))+r1*U(k,2:nx) ...
                -U(k-1,2:nx)+ht2*feval(f,xvec(2:nx),tvec(k));
end
%  plot the numerical solution
surf(xvec,tvec,U)
xlabel('x-axis')
ylabel('t-axis')
zlabel('The numerical solution')
s1=sprintf('h_t=%6.4f  h_x=%6.4f',ht,hx)
title(['Solution: ',s1])
