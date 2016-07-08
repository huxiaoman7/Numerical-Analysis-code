function value=polyeval(x,alpha,coeff,n);
%
% function value=polyeval(x,coeff,n)
%
% Evaluate a Taylor polynomial at the points given in x, with 
% alpha the point of expansion of the Taylor polynomial, and
% with n the degree of the polynomial.  The coefficients are to
% be given in coeff; and it is assumed there are n+1 entries in
% coeff with coeff(1) the constant term in the polynomial
%
value=coeff(n+1)*ones(size(x));
z=x-alpha;
for i=n:-1:1
  value = coeff(i) + z.*value;
end
