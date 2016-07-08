% Title: Plot a 3D graph of variation in a Taylor series approximation
% as the degree varies.  This program is for log(x) about x=1.

% This produces a surface plot to show the function values of a
% Taylor polynomial approximation for increasing values of the 
% degree of the approximation.  It also shows the values of the 
% errors for increasing degrees.

% Initialize
max_deg=10;                       % Set the maximum degree.
c=log_tay(max_deg);               % Call appropriate generator of Taylor coefficients.

% Do the Taylor series calculation over [a,b].
a=.1; b=2; step=.01;
x=a:step:b;                       % Generate values of independent variable.
L=length(x);

% Set up the z-values for the 3D-plot.
z=zeros(max_deg+1,L);             % Set up array z.
alpha=1;                          % Point of expansion.
for n=1:max_deg
  z(n,:)=polyeval(x,alpha,c,n);
end
z(max_deg+1,:)=log(x);            %Set true values of function.

% Set up xy-grid points
y=1:max_deg+1;
[X,Y]=meshgrid(x,y);

% Plot the function values.
mesh(X,Y,z),grid;
title('Graphs of Taylor polynomial for log(x)')    %Set title.
xlabel('The variable x')
ylabel('The degree n')
pause

% print
% print -deps log3d.eps

% Set up grid values and z-values for 3D-plotting of errors.
y=1:max_deg;
[X,Y]=meshgrid(x,y);
true=log(x);
err=zeros(max_deg,L);
for n=1:max_deg
  err(n,:)=z(n,:)-true;
end

% Plot the error values.
mesh(X,Y,err),grid;
title('Errors in Taylor polynomial for log(x)')   %Set title.
xlabel('The variable x')
ylabel('The degree n')

%print
%print -deps logerr3d.eps
