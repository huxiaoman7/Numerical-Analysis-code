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
