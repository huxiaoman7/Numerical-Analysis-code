function value=poly_odd(x,coeff,n);
%
% function value=poly_odd(x,coeff,n)
%
% Evaluate an "odd" Taylor polynomial at the points given in x,
% with n the degree of the polynomial.  The coefficients are to
% be given in coeff.  It is assumed that the numbers in coeff are
% are the coefficients of the odd-ordered terms in the polynomial.
%
m=ceil(n/2);
if n~=2*m-1
  'Error in poly_odd: the degree n must be an odd integer.'
  return
end
xsq = x.*x;
value=coeff(m)*ones(size(x));
for i=m-1:-1:1
  value = coeff(i) + xsq.*value;
end
value=x.*value;
