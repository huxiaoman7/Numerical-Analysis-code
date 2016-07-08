function coeff=sint_tay(n)
%
% function coeff = sint_tay(n)
%
% Evaluate the coefficients of the Taylor approximation of 
% degree n which approximates the "sine integral".  The input
% variable n must be an even integer.  The output vector coeff
% will have length m+1 where m=n/2.  This is because only the 
% even ordered coefficients are nonzero, and this must be
% considered in evaluating the associated Taylor polynomial.
%
m = double(int32(n/2));
if n ~= 2*m
  disp('Error in poly_even(x,coeff,n):')
  disp('The parameter n must be an even integer.')
  return
end
%
coeff = ones(m+1,1);
sign = 1;
fact = 1;
%
for i=2:m+1
  sign = -sign;
  d = 2*i-1;
  fact = fact*(d-1)*d;
  coeff(i) = sign/(fact*d);
end 
