function [nodes, fcn_values, div_diff_fcn] = chebyshev_interp(n)

% This creates an interpolant of order n to the function
% fcn(x) on [-1,1],  which is given below as a function 
% subprogram. The nodes are the Chebyshev zeroes of the 
% degree n+1 Chebyshev polynomial on [-1,1]. The program 
% gives two plots: first the true function and its 
% interpolant, and second, the error in the interpolation.

% Create the nodes and associated divided differences.
h = pi/(2*(n+1));
nodes = cos(h*[1:2:2*n+1]);
fcn_values = fcn(nodes);
div_diff_fcn = divdif(nodes,fcn_values);

% Create the points at which the functions are to be 
% graphed.
x_eval = -1:.002:1;
true_fcn = fcn(x_eval);
y_eval = interp(nodes,div_diff_fcn,x_eval);

% Create the graph of the function and its interpolant.
m = min([min(true_fcn),min(y_eval)]);
if m > 0
    m = .9*m;
else
    m = 1.1*m;
end
M = max([max(true_fcn),max(y_eval)]);
if M < 0
    M = .9*M;
else
    M = 1.1*M;
end
axis([-1.1,1.1,m,M])
hold on
plot(x_eval,true_fcn,'r','LineWidth',1)
plot(x_eval,y_eval,':')
legend('True function','Interpolant',0)
hold on
plot(nodes,fcn_values,'.','MarkerSize',6)
hold off

pause
clf

% Create the window for the graph of the error.
error = true_fcn - y_eval;
M = max(error);
if M < 0
    M = .9*M;
else
    M = 1.1*M;
end
m = min(error);
if m > 0
    m = .9*m;
else
    m = 1.1*m;
end
axis([-1.1,1.1,m,M])
hold on

% Create the graph of the error in the interpolant.
plot(x_eval,error,'r','LineWidth',1)
hold off

% Print the maximum error.
disp(['maximum error = ',num2str(max(abs(error)))])

function fval = fcn(x)

fval = exp(x);

function p_eval = interp(x_nodes,divdif_y,x_eval)
%
% This is a function
%         p_eval = interp(x_nodes,divdif_y,x_eval)
% It calculates the Newton divided difference form of 
% the interpolation polynomial of degree m-1, where the 
% nodes are given in x_nodes, m is the length of x_nodes, 
% and the divided differences are given in divdif_y.  The 
% points at which the interpolation is to be carried out 
% are given in x_eval; and on exit, p_eval contains the 
% corresponding values of the interpolation polynomial.
%
m = length(x_nodes);
p_eval = divdif_y(m)*ones(size(x_eval));
for i=m-1:-1:1
    p_eval = divdif_y(i) + (x_eval - x_nodes(i)).*p_eval;
end

function divdif_y = divdif(x_nodes,y_values)
%
% This is a function
%         divdif_y = divdif(x_nodes,y_values)
% It calculates the divided differences of the function 
% values given in the vector y_values, which are the values of
% some function f(x) at the nodes given in x_nodes.  On exit, 
%         divdif_y(i) = f[x_1,...,x_i],  i=1,...,m
% with m the length of x_nodes.  The input values x_nodes and
% y_values are not changed by this program.
%
divdif_y = y_values;
m = length(x_nodes);
for i=2:m
    for j=m:-1:i
        divdif_y(j) = (divdif_y(j)-divdif_y(j-1)) ...
                           /(x_nodes(j)-x_nodes(j-i+1));
    end
end
