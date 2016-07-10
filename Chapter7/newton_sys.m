function solution = newton_sys(x_init,err_tol,max_iterates)
%
% The calling sequence for newton_sys.m is
%    solution = newton_sys(x_init,err_tol,max_iterates)
% This solves a pair of two nonlinear equations in two unknowns,
%       f(x) = 0
% with f(x) a column vector of length 2. The definition of f(x)
% is to be given below in the function named fsys; and you 
% also need to give the Jacobian matrix for f(x) in the 
% function named deriv_fsys.
%
% x_init is a vector of length 2, and it is an initial guess 
% at the solution.
%
% The parameters err_tol and max_iterates are upper limits on
% the desired error in the solution and the maximum number of 
% iterates to be computed.
%
% Initialization.
x0=zeros(2,1);
for i=1:2
    x0(i) = x_init(i);
end

error = inf;
it_count = 0;

% Begin the main loop.
while error > err_tol & it_count < max_iterates
    it_count = it_count + 1;
    rhs = fsys(x0);
    A = deriv_fsys(x0);
    delta = A\rhs;
    x1 = x0 - delta;
    error = norm(delta,inf);
    % The following statement is an internal print to show
    % the course of the iteration.  It and the pause 
    % statement following it can be commented out.
    [it_count x1' error]
    pause
    x0 = x1;
end

% Return with the solution.
solution = x1;
if it_count == max_iterates
    disp('  ')
    disp('*** Your answers may possibly not satisfy your error tolerance.')
end

%%%%%%%%%%% Definition of functions %%%%%%%%%%%%%%%%%
function f_val = fsys(x)
%
% The equations being solved are
%       x(1)^2 + 4*x(2)^2 - 9 = 0
%    18*x(2) - 14*x(1)^2 + 45 = 0
%
f_val = [x(1)^2+4*x(2)^2-9, 18*x(2)-14*x(1)^2+45]';

function df_val = deriv_fsys(x)
%
% This defines the Jacobian matrix for the function
% given in fsys

df_val = [2*x(1), 8*x(2); -28*x(1), 18];
