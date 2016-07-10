function [lambda, eigenvec, diff, ratio] = ...
                         power_method(A,x_0,index,err_tol,max_iterates)
%
% The calling sequence for power_method.m is
%
% [lambda, eigenvec, diff, ratio] = power_method(A,x_0,index, ...
%                                                 err_tol,max_iterates)
%
% It implements the power method for finding the eigenvalue of A 
% which is largest in magnitude.  
%
% The vector x_0 is the user's initial guess at the eigenvector 
% corresponding to the eigenvalue of largest magnitude.  If the user has 
% no idea of how to choose x_0, then we suggest using
%      x_0 = rand(n,1)
% with n the order of A.  DO NOT use x_0 = 0.  
%
% The input parameter "index" is used to specify the method of 
% approximating the eigenvalue.  With index=0, we generate a random
% vector, call it v, and we use
%       lambda approximately (v'*w)/(v'*z)
% with Az=w.  For index in [1,n], we use
%       lambda approximately w(index)/z(index).
%
% The parameter err_tol is used in the error test.  However, note 
% this test needs to be made more sophisticated. It uses 
%       error approx lambda_present - lambda_previous
% at each step of the iteration, and this is not an accurate 
% estimate for slowly convergent iterations.
%
% The output variables are:
% lambda:    A vector of the power method iterations for the 
%            approximate eigenvalue of largest magnitude.
% eigenvec:  The corresponding normalized eigenvector.
% diff:      The differences of the successive values in lambda.
% ratio:     The ratios of the successive values in diff.
%            The values of diff and ratio are of use in examining
%            the rate of convergence of the power method.
%
% Initialization
n = size(A,1);				        % Find the order of A.
x_init = zeros(n,1);
lambda = zeros(max_iterates,1);
%
% Define x_init and require it to be a column vector.
for i=1:n
    x_init(i) = x_0(i);
end
%
[max_x,i_max] = max(abs(x_init));   % Find the max element in abs(x_init).
z = x_init/x_init(i_max);           % Normalize the initial vector to have
                                    %    an infinity norm of 1.
if index == 0                       % Generate special vector for 
   spec_vec = rand(1,n);            %    approximate eigenvalue computation.
end
%
format long
error = inf;
it_count = 0;
% Main loop.
while abs(error) > err_tol & it_count < max_iterates
    it_count = it_count + 1; 
    w = A*z;
    [max_w,i_max] = max(abs(w));
    if index == 0                   % Calculate approximate eigenvalue.
        lambda(it_count) = (spec_vec*w)/(spec_vec*z);
    else
        lambda(it_count) = w(index)/z(index); 
    end 
%    
    if it_count > 1
        error = lambda(it_count) - lambda(it_count-1);
    end
%    
    z = w/w(i_max);                 % Save approximate eigenvector.
end
%
if it_count == max_iterates
    disp('  ')
    disp('*** Your answers may possibly not satisfy your error tolerance.')
end
%
% Setting of output variables.
eigenvec = z;
lambda = lambda(1:it_count);
diff = zeros(it_count,1); 
ratio = zeros(it_count,1);
diff(2:it_count) = lambda(2:it_count)-lambda(1:it_count-1);
ratio(3:it_count) = diff(3:it_count)./diff(2:it_count-1);
