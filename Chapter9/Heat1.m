function U = Heat1(a,f,u0,g1,g2,L,T,nx,nt)
%
% function U = Heat1(a,f,u0,g1,g2,L,T,nx,nt)
%
% The forward difference scheme for solving the initial-boundary value
% problem of the heat equation u_t = a u_{xx} + f for x in [0,L], and
% t in [0,T].
% Input
%      a: the coefficient of the u_{xx} term
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
r = a*ht/(hx*hx)
r1 = 1-2*r;
% stability test
if r > 0.5
  disp('The ratio is too big; the scheme is unstable!')
end
% generate grid point vectors
xvec = hx*(0:nx);
tvec = ht*(0:nt);
% initialization
U = zeros(nt1,nx1);
U(1,:) = feval(u0,xvec(:)');
% time advancing
U(2:nt1,1) = feval(g1,(1:nt)*ht);
U(2:nt1,nx1) = feval(g2,(1:nt)*ht);
for k = 1:nt
  U(k+1,2:nx) = r*(U(k,1:nx-1)+U(k,3:nx+1))+r1*U(k,2:nx) ...
                +ht*feval(f,xvec(2:nx),tvec(k));
end
% plot the numerical solution
surf(xvec,tvec,U)
xlabel('x-axis')
ylabel('t-axis')
zlabel('The numerical solution')
s1 = sprintf('h_t=%6.4f  h_x=%6.4f',ht,hx)
title(['Solution: ',s1])
