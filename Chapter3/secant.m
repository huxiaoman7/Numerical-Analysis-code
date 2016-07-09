
function root = secant(x0,x1,error_bd,max_iterate,index_f)
%
% function secant(x0,x1,error_bd,max_iterate,index_f)
%
% This implements the secant method for solving an 
% equation f(x) = 0. The function f(x) is given below.
%
% The parameter error_bd is used in the error test for the 
% accuracy of each iterate.  The parameter max_iterate is an 
% upper limit on the number of iterates to be computed. Two 
% initial guesses, x0 and x1, must also be given.
%
% For the given function f(x), an example of a calling sequence 
% might be the following:
%    root = secant(x0,x1,1.0E-12,10,1)
% The parameter index_f specifies the function to be used.
%
% The program prints the iteration values
%      iterate_number, x, f(x), error
% The value of x is the most current initial guess, called
% previous_iterate here, and it is updated with each iteration. 
% The value of error is 
%    error = newly_computed_iterate - previous_iterate
% and it is an estimated error for previous_iterate.
% Tap the carriage return to continue with the iteration. 

format short e
error = 1;
fx0 = f(x0,index_f);
it_count = 0;
iteration = [it_count x0 fx0]
while abs(error) > error_bd & it_count <= max_iterate
    it_count = it_count + 1;
    fx1 = f(x1,index_f);
    if fx1 - fx0 == 0
        disp('f(x1) = f(x0); Division by zero;  Stop')
        return
    end
    x2 = x1 - fx1*(x1-x0)/(fx1-fx0);
    error = x2 - x1;
%   Internal print of secant method. Tap the carriage
%   return key to continue the computation.
    iteration = [it_count x1 fx1 error]
    pause
    x0 = x1;
    x1 = x2;
    fx0 = fx1;
end

if it_count > max_iterate
    disp('The number of iterates calculated exceeded') 
    disp('max_iterate.  An accurate root was not')
    disp('calculated.')
else
    format long
    root = x2
    format short e
    error
    format short
    it_count
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function value = f(x,index)

% function to define equation for rootfinding problem.

switch index
case 1        
    value = x.^6 - x - 1;
case 2
    value = x - exp(-x);
end

