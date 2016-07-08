function coeff = log_tay(n);
%
% function coeff = log_tay(n)
%
% Produce the Taylor coefficients for the function log(x) 
% when expanded about the point a=1.  The coefficients will
% be returned in the array coeff, which will have length n+1.
%
coeff = zeros(n+1,1);
sign = 1;
for i=1:n
  coeff(i+1) = sign/i;
  sign = -sign;
end
