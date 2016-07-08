function ans = plot_sint_total

% TITLE: Plot Taylor polynomials for the "sine integral" 
%        about x = 0.
%
% This plots several Taylor polynomials and their errors
% for increasing degrees.  The particular function being
% approximated is Sint(x) on [0,b], with x = 0 the point of
% expansion for creating the Taylor polynomials.  We plot the 
% Taylor polynomials for several degrees, which can be altered
% by the user.  Note that the function being plotted is symmetric
% about x=0.
%
% The Taylor polynomials in this case contain only terms of even 
% degree.  Such polynomials are called "even polynomials", and 
% they are symmetric about x=0.
%
% TO THE STUDENT: run this program for various input values of "b".
% Also, change the values given in the vector "degree" given below, to
% experiment with different degrees of Taylor polynomial approximations.
%
% Initialize
b = input('Give the number b defining the interval [0,b] ');
h = b/200;
x = 0:h:b;
max_degree = 20;
%
% Produce the Taylor coefficients for the "sine integral" function "sint".
c = sint_tay(max_degree);
%
% Specify the four values of degree to be considered.  They must all be even, 
% and they must be less than or equal to max_degree.
degree=[2,4,6,8];
if max(degree) > max_degree
    fprintf('Some value of degree is greater than max_degree = %2.0f\n',...
             max_degree)
    return
end
%
% Initialize the array to contain the polynomial values.  Row #i is to 
% contain the values for the polynomial of degree=degree(i).
p = zeros(4,length(x));
%
% Calculate the Taylor polynomials
for i=1:4
    p(i,:) = poly_even(x,c,degree(i));
end
%
% Initialize for plotting the Taylor polynomials
hold off
clf
axis([0,b,0,1])
hold on
%
% Plot the Taylor polynomials
plot(x,p(1,:),x,p(2,:),':',x,p(3,:),'--',x,p(4,:),'-.')
plot([0,b],[0,0])
plot([0 0],[0 1])
title('Taylor approximations of Sint(x)')
text(1.025*b,0,'x')
text(0,1.03,'y')
legend(strcat('degree = ',int2str(degree(1))),...
       strcat('degree = ',int2str(degree(2))),...
       strcat('degree = ',int2str(degree(3))),...
       strcat('degree = ',int2str(degree(4))),0)

function coeff=sint_tay(n)
%
% function coeff = sint_tay(n)
%
% Evaluate the coefficients of the Taylor approximation of 
% degree n which approximates the "sine integral".  The input
% variable n must be an even integer.  The output vector coeff
% will have length m+1 where m=n/2.  This is because only the 
% even ordered coefficients are nonzero, and this must be
% considered in evaluating the associated Taylor polynomial.
%
m = double(int32(n/2));
if n ~= 2*m
  disp('Error in poly_even(x,coeff,n):')
  disp('The parameter n must be an even integer.')
end
%
coeff = ones(m+1,1);
sign = 1;
fact = 1;
%
for i=2:m+1
  sign = -sign;
  d = 2*i-1;
  fact = fact*(d-1)*d;
  coeff(i) = sign/(fact*d);
end 

function value = poly_even(x,coeff,n);
%
% function value = poly_even(x,coeff,n)
%
% Evaluate an "even" Taylor polynomial at the points 
% given in x, with n the degree of the polynomial.  
% The coefficients are to be given in coeff.  It is 
% assumed that the numbers in coeff are are the non-zero  
% coefficients of the even-ordered terms in the polynomial.  
% The input parameter n must be an even integer.
%
m = double(int32(n/2));
if n ~= 2*m
  disp('Error in poly_even(x,coeff,n):')
  disp('The parameter n must be an even integer.')
end
%
xsq = x.*x;
value = coeff(m+1)*ones(size(x));
%
for i = m:-1:1
  value = coeff(i) + xsq.*value;
end
