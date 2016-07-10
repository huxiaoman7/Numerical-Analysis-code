function [x,lu,piv] = GEpivot(A,b)
%
% function [x,lu,piv] = GEpivot(A,b)
%
% This program employs the Gaussian elimination method with partial
% pivoting to solve the linear system Ax=b.  
%
% Input
%    A: coefficient square matrix
%    b: right side column vector
%
% Output
%    x: solution vector
%    lu: a matrix whose upper triangular part is the upper 
%        triangular matrix resulting from the Gaussian elimination
%        with partial pivoting, and whose strictly lower triangular
%        part stores the multipliers.  In other words, it stores
%        the LU factorization of the matrix A with row permutations
%        determined by the pivoting vector piv.  
%    piv: a pivoting vector recording the row permutations.

% check the order of the matrix and the size of the vector
[m,n] = size(A);
if m ~= n
   error('The matrix is not square.')
end
m = length(b);
if m ~= n
   error('The matrix and the vector do not match in size.')
end

% initialization of the pivoting vector
piv = (1:n)';

% elimination step
for k = 1:n-1
   % Find the maximal element in the pivot column, below the 
   % pivot position, along with the index of that maximal element.
   
   [col_max  index] = max(abs(A(k:n,k)));
   index = index + k-1;
   
   if index ~= k
      % Switch rows k and index, in columns k thru n.  Do similarly
      % for the right-hand side b.
      
      tempA = A(k,k:n);
      A(k,k:n) = A(index,k:n);
      A(index,k:n) = tempA;
      
      tempb = b(k);
      b(k) = b(index);
      b(index) = tempb;
      
      temp = piv(k);
      piv(k) = piv(index);
      piv(index) = temp;
   end
   
   % Form the needed multipliers and store them into the pivot
   % column, below the diagonal.
   A(k+1:n,k) = A(k+1:n,k)/A(k,k);
   
   % Carry out the elimination step, first modifying the matrix,
   % and then modifying the right-hand side.
   for i = k+1:n
      A(i,k+1:n) = A(i,k+1:n) - A(i,k)*A(k,k+1:n);
   end
   b(k+1:n) = b(k+1:n) - A(k+1:n,k)*b(k);
end

% Solve the upper triangular linear system.
x = zeros(n,1);
x(n) = b(n)/A(n,n);
for i = n-1:-1:1
   x(i)=(b(i)-A(i,i+1:n)*x(i+1:n))/A(i,i);
end

% Record the LU factorization with row permutation.
lu = A;
