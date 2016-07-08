function value = sint(x)
%
%   This evaluates an approximation to
%
%                 1     x    1
%       Sint(X) = - INTEGRAL - Sin(t)*dt
%                 x     0    t
%
%   For -1 .LE. x .LE. 1.  the approximation is the
%   degree 8 Taylor polynomial for Sint(x).  Its
%   maximum error is 5.0E-9, provided the computer
%   arithmetic allows an error of that size.
%

coeff = [1, -1/18, 1/600, -1/35280, 1/3265920];
degree = length(coeff)-1;

u = x.*x;

value = coeff(degree+1)*ones(size(x));

for j=degree:-1:1;
    value = coeff(j) + u.*value;
end
    
