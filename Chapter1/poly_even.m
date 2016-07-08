function value = poly_even(x,coeff,n);
%
% function value = poly_even(x,coeff,n)
%
% Evaluate an "even" Taylor polynomial at the points 
% given in x, with n the degree of the polynomial.  
% The coefficients are to be given in coeff.  It is 
% assumed that the numbers in coeff are are the non-zero  
% coefficients of the even-ordered terms in the polynomial.  
% The input parameter n must be an even integer.
%
m = ceil(n/2);
if n ~= 2*m
  disp('Error in poly_even(x,coeff,n):')
  disp('The parameter n must be an even integer.')
  return
end
%
xsq = x.*x;
value = coeff(m+1)*ones(size(x));
%
for i = m:-1:1
  value = coeff(i) + xsq.*value;
end
