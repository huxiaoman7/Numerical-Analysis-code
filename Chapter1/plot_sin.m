% TITLE: Plot Taylor polynomials for sin(x) about x = 0
%
% This plots several Taylor polynomials and their errors
% for increasing degrees.  The particular function being
% approximated is sin(x) on [0,b], with x = 0 the point of
% expansion for creating the Taylor polynomials.  First  
% we plot the Taylor polynomials and then we plot the errors.
%
% Note that the Taylor polynomials in this case contain only
% terms of odd degree.  Such polynomials are called "odd
% polynomials".
%
% Initialize
b = input('Give the number b defining the interval [0,b] ');
h = b/200;
x = 0:h:b;
max_deg = 7;

% Produce the Taylor coefficients for the function sin(x) 
c = sin_tay(max_deg);

% Calculate the Taylor polynomials
p1 = poly_odd(x,c,1);
p3 = poly_odd(x,c,3);
p5 = poly_odd(x,c,5);
p7 = poly_odd(x,c,7);

% Calculate the errors in the Taylor polynomials
true = sin(x);
err1 = true-p1;
err3 = true-p3;
err5 = true-p5;
err7 = true-p7;

% Initialize for plotting the Taylor polynomials
hold off
clf
m = 1.1*min([min(p1),min(p3),min(p5),min(p7)]);
M = 1.1*max([max(p1),max(p3),max(p5),max(p7)]);
axis([0,b,m,M])
hold on

% Plot the Taylor polynomials
plot(x,true,x,p1,':',x,p3,'--',x,p5,'-.');
title('Taylor approximations of sin(x)')
plot([0,b],[0,0])
plot([0 0 ],[0 M/1.2])
axis('equal')
legend('sin(x)','degree = 1','degree = 3','degree = 5')
pause

% print -deps sin2d.eps
% print

% Initialize for plotting the errors
hold off
clf
m = 1.2*min([min(err1),min(err3),min(err5),min(err7)]);
M = 1.2*max([max(err1),max(err3),max(err5),max(err7)]);
axis([0,b,m,M])
hold on

% Plot the errors
plot(x,err1,x,err3,':',x,err5,'--',x,err7,'-.');
title('Errors in Taylor approximations of sin(x)')
legend('degree = 1','degree = 3','degree = 5','degree = 7')
hold off

% print -deps sinerr2d.eps
% pause
% print
