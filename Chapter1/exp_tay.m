function coeff=exp_tay(n);
%
% function coeff=exp_tay(n)
%
% Produce the Taylor coefficients for the function exp(x) 
% when expanded about the point a=0.  The coefficients will
% be returned in the array coeff, which will have length n+1.
%
coeff=ones(n+1,1);
fact=1;
for i=1:n
  fact=i*fact;
  coeff(i+1)=1/fact;
end
