function [x, iflag, itnum] = Jacobi(A,b,x0,delta,max_it)
%
% function [x, iflag, itnum] = Jacobi(A,b,x0,delta,max_it)
%
% A program implementing the Jacobi iteration method to solve
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
% create a vector with diagonal elements of A
diagA = diag(A);
% modify A to make its diagonal elements zero
A = A-diag(diag(A));
% iteration
while k < max_it
  k = k+1;
  x = (b-A*x0)./diagA
  relerr = norm(x-x0,inf)/(norm(x,inf)+eps);
  x0 = x;
  if relerr < delta
    break
  end
end
%
itnum = k;
if (itnum == max_it)
  iflag = -1
end
