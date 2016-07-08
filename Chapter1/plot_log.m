% TITLE: Plot Taylor polynomials for log(x) about x = 1
%
% This plots several Taylor polynomials and their errors
% for increasing degrees.  The particular function being
% approximated is log x and it is being expanded about 
% x = 1. The interval of approximation is [a,b] with   
% 0 < a < 1 and b > 1.  First we plot the Taylor polynomials and   
% then we plot the errors.
%
% Initialize
a = input('For the interval [a,b], give a with 0 < a < 1. ');
b = input('Give b > 1 ');
h = (b-a)/200;
x = a:h:b;
max_deg = 4;
c = zeros(max_deg+1,1);
sign = 1;
for i = 1:max_deg
  c(i+1) = sign/i;
  sign = -sign;
end

% Calculate the Taylor polynomials
p1 = polyeval(x,1,c,1);
p2 = polyeval(x,1,c,2);
p3 = polyeval(x,1,c,3);
p4 = polyeval(x,1,c,4);

% Calculate the errors in the Taylor polynomials
true = log(x);
err1 = true-p1;
err2 = true-p2;
err3 = true-p3;
err4 = true-p4;

% Initialize for plotting the Taylor polynomials
hold off
clf
m = 1.1*min([min(p1),min(p2),min(p3),min(p4)]);
M = 1.1*max([max(p1),max(p2),max(p3),max(p4)]);
axis([a,b,m,M])
hold on

% Plot the Taylor polynomials
plot(x,true,x,p2,':',x,p3,'--',x,p4,'-.');
plot([0 b],[0 0])
plot([0 0],[m/1.1 M/1.1])
title('Taylor approximations of log(x)')
text(b+3*h,0,'x')
text(0,M,'y')
axis('equal')
legend('log(x)','degree = 2','degree = 3','degree = 4')
pause

% print -deps log2d.eps
% print

% Initialize for plotting the errors
hold off
clf
m = 1.1*min([min(err1),min(err2),min(err3),min(err4)]);
M = 1.1*max([max(err1),max(err2),max(err3),max(err4)]);
axis([0,b,m,M])
hold on

% Plot the errors
plot(x,err1,x,err2,':',x,err3,'--',x,err4,'-.');
plot([0 b],[0 0])
plot([0 0],[m/1.2 M/1.2])
title('Errors in Taylor approximations of log(x)')
text(b+3*h,0,'x')
text(0,M+.05*max(M,abs(m)),'y')
legend('degree = 1','degree = 2','degree = 3','degree = 4')

% pause
% print -deps logerr2d.eps
% print
