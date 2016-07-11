function [x,y] = AB2(x0,y0,x_end,h,fcn)
%
% function [x,y]=AB2(x0,y0,x_end,h,fcn)
%
% Solve the initial value problem
%    y' = f(x,y),  x0 <= x <= b,  y(x0)=y0
% Use Adams-Bashforth formula of order 2 with a stepsize of h.  
% Euler's method is used for the value y1.  The user must 
% supply a program with some name, say deriv, and a first 
% line of the form
%   function ans=deriv(x,y)
% A sample call would be
%   [t,z]=AB2(t0,z0,b,delta,'deriv')
%
% Output:
% The routine AB2 will return two vectors, x and y.
% The vector x will contain the node points
%   x(1)=x0, x(j)=x0+(j-1)*h, j=1,2,...,N
% with 
%   x(N) <= x_end-h,  x(N)+h > x_end-h
% The vector y will contain the estimates of the solution Y 
% at the node points in x.
%
n = fix((x_end-x0)/h)+1;
x = linspace(x0,x0+(n-1)*h,n)';
y = zeros(n,1);
y(1) = y0;
ft1 = feval(fcn,x(1),y(1));
y(2) = y(1)+h*ft1;
for i = 3:n
  ft2 = feval(fcn,x(i-1),y(i-1));
  y(i) = y(i-1)+h*(3*ft2-ft1)/2;
  ft1 = ft2;
end



