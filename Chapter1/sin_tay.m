function coeff=sin_tay(n)
%
% function coeff=sin_tay(n)
%
% Evaluate the coefficients of the Taylor approximation of 
% degree n which approximates the sine function.  The input
% variable n must be an odd integer.  The output vector coeff
% will have length m where n=2m-1.  This is because only the 
% odd ordered coefficients are nonzero, and this must be
% considered in evaluating the associated Taylor polynomial.
%
m = ceil(n/2);
if n ~= 2*m-1 
  'Error in sin_tay(n): the variable n must be an odd integer.'
end

coeff=ones(m,1);
sign = 1;
fact = 1;

for i=2:m
  sign = -sign;
  fact = fact*(2*i-2)*(2*i-1);
  coeff(i) = sign/fact;
end 
