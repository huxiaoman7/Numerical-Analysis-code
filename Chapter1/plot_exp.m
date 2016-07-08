% TITLE: Plot Taylor polynomials for exp(x) about x = 0
%
% This plots several Taylor polynomials and their errors
% for increasing degrees.  The particular function being
% approximated is exp(x) on [-b,b].  First we plot the 
% Taylor polynomials and then we plot the errors.
%
% Initialize
b = input('Give the number b defining the interval [-b,b] ');
h = b/100;
x = -b:h:b;
max_deg = 4;

% Does the user want equal scaling?
disp('The phrase "equal scaling" means the horizontal')
disp('and vertical axes have the same unit length.')
disp('Do you want equal scaling?')
query = input('Give 0 for "no" and 1 for "yes". ');

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

% Initialize for plotting the Taylor polynomials
hold off
clf
m = min([min(p1),min(p2),min(p3),min(p4),min(true)]);
if m > 0
    m = .9*m;
else
    m = 1.1*m;
end
M = max([max(p1),max(p2),max(p3),max(p4),max(true)]);
if M < 0
    M = .9*M;
else
    M = 1.1*M;
end
axis([-b,b,m,M])
hold on

% Plot the Taylor polynomials
plot(x,true,'k',x,p1,'b:',x,p2,'g--',x,p3,'r-.');
title('Taylor approximations of e^x')
plot([-b,b],[0,0])
plot([0 0 ],[m M/1.2])
if query == 1
    axis('equal')
end
legend('exp(x)','degree = 1','degree = 2','degree = 3',0)
pause

% print -deps exp2d.eps
% print

% Initialize for plotting the errors
hold off
clf
m = min([min(err1),min(err2),min(err3),min(err4)]);
if m > 0
    m = .9*m;
else
    m = 1.1*m;
end
M = max([max(err1),max(err2),max(err3),max(err4)]);
if M < 0
    M = .9*M;
else
    M = 1.1*M;
end
axis([-b,b,m,M])
hold on

% Plot the errors
plot(x,err1,'b',x,err2,'g:',x,err3,'r--',x,err4,'c-.');
title('Errors in Taylor approximations of e^x')
legend('degree = 1','degree = 2','degree = 3','degree = 4',0)
hold off

% print -deps experr2d.eps
% pause
% print
