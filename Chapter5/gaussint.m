function [val,bp,wf]=gaussint(a,b,n,index_f)

% [val,bp,wf]=gaussint(fun,a,b,n) integrates
% a function from a to b using an n-point
% Gauss rule which is exact for a polynomial
% of degree 2*n-1. Concepts on page 93 of
% 'Methods of Numerical Integration' by
% Philip Davis and Philip Rabinowitz yield
% the base points and weight factors.
%
% a,b   -  integration limits
% n     -  order of formula (default is 20)
% val   -  value of the integral
% bp,wf -  Gauss base points and weight factors on [-1,1]
% index_f - index of function to be integrated.
 
if nargin < 4, n=20; end
u=(1:n-1)./sqrt((2*(1:n-1)).^2-1);
[vc,bp]=eig(diag(u,-1)+diag(u,1));
[bp,k]=sort(diag(bp));
wf=2*vc(1,k)'.^2; 
x=(a+b)/2+((b-a)/2)*bp;
f=fcn(x,index_f)*(b-a)/2; 
val=wf(:)'*f(:);


function f_value = fcn(x,index)
%
% This defines the integrand.

switch index
case 1
    f_value = exp(-x.^2);
case 2
    f_value = 1 ./(1+x.^2);
case 3
    f_value = 1 ./(2+sin(x));
case 4
    f_value = exp(cos(x));
end
