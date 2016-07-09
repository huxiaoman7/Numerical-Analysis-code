function root=bisect(a0,b0,ep,max_iterate,index_f)
%
% function bisect(a0,b0,ep,max_iterate,index_f)
%
% This is the bisection method for solving an equation f(x)=0.
%
% The function f is defined below by the user. The function f is 
% to be continuous on the interval [a0,b0], and it is to be of 
% opposite signs at a0 and b0.  The quantity ep is the error 
% tolerance.  The routine guarantees this as an error bound 
% provided: (1) the restrictions on the initial interval are 
% correct, and (2) ep is not too small when the machine epsilon 
% is taken into account.  Most of these conditions are not 
% checked in the program!  The parameter max_iterate is an upper 
% limit on the number of iterates to be computed.
%
% For the given function f(x), an example of a calling sequence 
% might be the following:
%    root = bisect(1,1.5,1.0E-6,10,1)
% The parameter index_f specifies the function to be used.
%
% The following will print out for each iteration the values of
%      count, a, b, c, f(c), (b-a)/2
% with c the current iterate and (b-a)/2 the error bound for c.
% The variable count is the index of the current interate.  Tap 
% the carriage return to continue with the iteration. 

if a0 >= b0
    disp('a0 < b0 is not true.  Stop!')
    return
end

format short e
a = a0; b = b0;
fa = f(a,index_f); fb = f(b,index_f);

if sign(fa)*sign(fb) > 0
    disp('f(a0) and f(b0) are of the same sign.  Stop!')
    return
end

c = (a+b)/2;
it_count = 0;
while b-c > ep & it_count < max_iterate
    it_count = it_count + 1;
    fc = f(c,index_f);
%   Internal print of bisection method. Tap the carriage
%   return key to continue the computation.
    iteration = [it_count a b c fc b-c]
    if sign(fb)*sign(fc) <= 0
        a = c;
        fa = fc;
    else
        b = c;
        fb = fc;
    end
    c = (a+b)/2;
    pause
end

format long
root = c
format short e
error_bound = b-c
format short
it_count

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function value = f(x,index)

% function to define equation for rootfinding problem.

switch index
case 1        
    value = x.^6 - x - 1;
case 2
    value = x - exp(-x);
end

