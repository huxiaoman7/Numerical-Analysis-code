function [x, iflag, itnum] = GS(A,b,x0,delta,max_it)
%
% function [x, iflag, itnum] = GS(A,b,x0,delta,max_it)
%
% A program implementing the Gauss-Seidel iteration method to solve
% the linear system Ax=b.
%
% Input
%   A: square coefficient matrix
%   b: right side vector
%   x0: initial guess 
%   delta: error tolerance for the relative difference between 
%          two consecutive iterates
%   max_it: maximum number of iterations to be allowed
%
% Output
%   x: numerical solution vector
%   iflag: 1 if a numerical solution satisfying the error 
%            tolerance is found within max_it iterations
%          -1 if the program fails to produce a numerical 
%             solution in max_it iterations
%   itnum: the number of iterations used to compute x
%

% initialization
n = length(b);
iflag = 1;
k = 0;
x = x0;
% iteration
while k < max_it
  k = k+1;
  x(1) = (b(1)-A(1,2:n)*x0(2:n))/A(1,1);
  for i = 2:n
    if i < n
      x(i) = (b(i)-A(i,1:i-1)*x(1:i-1) ...
              -A(i,i+1:n)*x0(i+1:n))/A(i,i);
    else
      x(n) = (b(n)-A(n,1:n-1)*x(1:n-1))/A(n,n);
    end
  end
  relerr = norm(x-x0,inf)/(norm(x,inf)+eps);
  x0 = x
  if relerr < delta
    break
  end
end
%
itnum = k;
if (itnum == max_it)
  iflag = -1
end
