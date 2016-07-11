function [x,y] = eulersys(x0,y0,x_end,h,fcn)
%
% function [x,y]=eulersys(x0,y0,x_end,h,fcn)
%
% Solve the initial value problem of a system
% of first order equations
%    y' = f(x,y),  x0 <= x <= b,  y(x0)=y0
% Use Euler's method with a stepsize of h.  
% The user must supply a program to compute the
% right hand side function with some name, say
% deriv, and a first line of the form
%   function ans=deriv(x,y)
% A sample call would be
%   [t,z]=eulersys(t0,z0,b,delta,'deriv')
%
% The program automatically determines the
% number of equations from the dimension of
% the initial value vector y0.
%
% Output:
% The routine eulersys will return a vector x
% and matrix y. The vector x will contain the
% node points in [x0,x_end]:
%   x(1)=x0, x(j)=x0+(j-1)*h, j=1,2,...,N
% The matrix y is of size N by m, with m the
% number of equations. The i-th row y(i,:) will
% contain the estimates of the solution Y 
% at the node points in x(i).
%
m = length(y0);
n = fix((x_end-x0)/h)+1;
x = linspace(x0,x0+(n-1)*h,n)';
y = zeros(n,m);
y(1,:) = y0;
for i = 2:n
  y(i,:) = y(i-1,:) + h*feval(fcn,x(i-1),y(i-1,:));
end
