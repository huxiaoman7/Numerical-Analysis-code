function [x,y]=trapezoidal(x0,y0,x_end,h,fcn,tol)
%
% function [x,y]=trapezoidal(x0,y0,x_end,h,fcn,tol)
%
% Solve the initial value problem
%    y' = f(x,y),  x0 <= x <= b,  y(x0)=y0
% Use trapezoidal method with a stepsize of h.  The user must 
% supply an m-file to define the derivative f, with some name, 
% say 'deriv.m', and a first line of the form
%   function ans=deriv(x,y)
% tol is the user supplied bound on the difference between successive
% values of the trapezoidal iteration.
% A sample call would be
%   [t,z]=trapezoidal(t0,z0,b,delta,'deriv',1e-3)
%
% Output:
% The routine trapezoidal will return two vectors, x and y.
% The vector x will contain the node points
%   x(1)=x0, x(j)=x0+(j-1)*h, j=1,2,...,N
% with 
%   x(N) <= x_end,  x(N)+h > x_end
% The vector y will contain the estimates of the solution Y 
% at the node points in x.
%

% Initialize.
n=fix((x_end-x0)/h)+1;
x=linspace(x0,x0+(n-1)*h,n)';
y=zeros(n,1);
y(1)=y0;
i=2;
% advancing
while i<= n
  fyt=feval(fcn,x(i-1),y(i-1));
%
% Euler estimate
%
  yt1=y(i-1)+h*fyt;
% trapezoidal iteration
  count=0;
  diff=1;
  while diff > tol & count < 10
    yt2=y(i-1) + h*(fyt+feval(fcn,x(i),yt1))/2;
    diff = abs(yt2-yt1);
    yt1=yt2;
    count = count +1;
  end
  if count >= 10
    disp('Not converging after 10 steps at x = ')
    fprintf('%5.2f\n', x(i))
  end
  y(i) = yt2;
  i=i+1;
end
