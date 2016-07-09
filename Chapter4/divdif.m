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
