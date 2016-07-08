% TITLE: Evaluate Taylor polynomials for exp(x) about x = 0
%
% This evaluates several Taylor polynomials and their errors
% for increasing degrees.  The particular function being
% approximated is exp(x) on [-b,b].  

% Initialize
b = input('Give the number b defining the interval [-b,b] ');
h = b/10;
x = -b:h:b;
max_deg = 4;

% Produce the Taylor coefficients for the function exp(x) when
% expanded about the point a = 0.  The coefficients are stored 
% in the array c, which will have length max_deg+1.
c = ones(max_deg+1,1);
fact = 1;
for i = 1:max_deg
  fact = i*fact;
  c(i+1) = 1/fact;
end

% Calculate the Taylor polynomials
p1 = polyeval(x,0,c,1);
p2 = polyeval(x,0,c,2);
p3 = polyeval(x,0,c,3);
p4 = polyeval(x,0,c,4);

% Calculate the errors in the Taylor polynomials
true = exp(x);
err1 = true-p1;
err2 = true-p2;
err3 = true-p3;
err4 = true-p4;

diary exp_taylor
disp('    x       exp(x)      err1          err2          err3          err4')
[x',true',err1',err2',err3',err4']
diary off
