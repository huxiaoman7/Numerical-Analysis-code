function [integral,difference,ratio]=trapezoidal(a,b,n0,index_f)
%
%  function [integral,difference,ratio]=trapezoidal(a,b,n0,index_f)
%
%  This uses the trapezoidal rule with n subdivisions to  
%  integrate the function f over the interval [a,b].  The 
%  values of n used are
%        n = n0,2*n0,4*n0,...,256*n0
%  The corresponding numerical integrals are returned in the 
%  vector integral.  The differences of successive numerical 
%  integrals are returned in the vector difference:
%        difference(i) = integral(i)-integral(i-1),  i=2,...,9
%  The entries in ratio give the rate of decrease in these
%  differences.
%
%  In using this program, define the integrand using the 
%  function given below.  The parameter index_f allows the user 
%  to do calculations with multiple integrands.

%  Initialize output vectors.
integral = zeros(9,1);
difference = zeros(9,1);
ratio = zeros(9,1);

% Initialize for trapezoidal rule.
sumend = (f(a,index_f) +f(b,index_f))/2;
sum = 0;

%  Initialize for case of n0 > 2.
if(n0 > 2)
    h = (b-a)/n0;
    for i=2:2:n0-2
        sum = sum + f(a+i*h,index_f);
    end
end

%  Calculate the numerical integrals, doing each
%  by appropriately modifying the preceding case.
for i=1:9
    n = n0*2^(i-1);
    h = (b-a)/n;
    for k=1:2:n-1
        sum = sum + f(a+k*h,index_f);
    end
    integral(i) = h*(sumend + sum);
end

%  Calculate the differences of the successive 
%  trapezoidal rule integrals and the ratio
%  of decrease in these differences.
difference(2:9)=integral(2:9)-integral(1:8);
ratio(3:9)=difference(2:8)./difference(3:9);

function f_value = f(x,index)
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
