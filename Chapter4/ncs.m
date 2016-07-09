function y_eval = ncs(x_nodes,y_nodes,x_eval)

m=length(x_nodes);
n=m-2;
a=zeros(1,n); b=zeros(1,n); c=zeros(1,n); f=zeros(1,n);
for i=1:n
    b(i)=(x_nodes(i+2)-x_nodes(i))/3;
end
for i=2:n
    a(i)=(x_nodes(i+1)-x_nodes(i))/6;
end
for i=1:n-1
    c(i)=(x_nodes(i+2)-x_nodes(i+1))/6;
end
for i=1:n
    f(i)=(y_nodes(i+2)-y_nodes(i+1))/(x_nodes(i+2)-x_nodes(i+1)) ...
            -(y_nodes(i+1)-y_nodes(i))/(x_nodes(i+1)-x_nodes(i));
end

[soln, alpha, beta, ier] = tridiag(a,b,c,f,n,0);

M=zeros(m);
M(2:m-1)=soln;

len_x=length(x_eval);
for i=1:len_x
    x=x_eval(i);
    for j=2:m
        if x_nodes(j-1) <= x & x <= x_nodes(j) 
            y_eval(i)=(M(j-1)*(x_nodes(j)-x)^3+M(j)*(x-x_nodes(j-1))^3)/(6*(x_nodes(j)-x_nodes(j-1))) ...
                +(y_nodes(j-1)*(x_nodes(j)-x)+y_nodes(j)*(x-x_nodes(j-1)))/(x_nodes(j)-x_nodes(j-1)) ...
                -(x_nodes(j)-x_nodes(j-1))*(M(j-1)*(x_nodes(j)-x)+M(j)*(x-x_nodes(j-1)))/6;
        end
    end
end

function [x, alpha, beta, ier] = tridiag(a,b,c,f,n,iflag)

% function [x, alpha, beta, ier] = tridiag(a,b,c,f,n,iflag)
%
% Solve a tridiagonal linear system M*x=f
%
% INPUT:
% The order of the linear system is given in n.
% The subdiagonal, diagonal, and superdiagonal of M are given
% by the arrays a,b,c, respectively.  More precisely,
%     M(i,i-1) = a(i), i=2,...,n
%     M(i,i)   = b(i), i=1,...,n
%     M(i,i+1) = c(i), i=1,...,n-1
% iflag=0 means that the original matrix M is given as specified above.
% iflag=1 means that the LU factorization of M is already known and is
% stored in a,b,c.  This will have been accomplished by a previous call
% to this routine.  In that case, the vectors alpha and beta should 
% have been substituted for a and b in the calling sequence.
%
% OUTPUT:
% Upon exit, the LU factorization of M is already known and is stored
% in alpha,beta,c.  The solution x is given as well.
% ier=0 means the program was completed satisfactorily.
% ier=1 means that a zero pivot element was encountered and the 
% solution process was abandoned.

a(1) = 0;
if iflag == 0
   % Compute LU factorization of matrix M.
   for j=2:n
      if b(j-1) == 0
         ier = 1;
         return
      end
      a(j) = a(j)/b(j-1);
      b(j) = b(j) - a(j)*c(j-1);
   end
   if b(n) == 0 
      ier = 1;
      return
   end
end

% Compute solution x to M*x = f.
% Do forward substitution to solve lower triangular system.
for j=2:n
   f(j) = f(j) - a(j)*f(j-1);
end

% Do backward substitution to solve upper triangular system.
f(n) = f(n)/b(n);
for j=n-1:-1:1
   f(j) = (f(j) - c(j)*f(j+1))/b(j);
end

% Set output variables.
ier = 0;
x = f;
alpha = a; beta = b;
