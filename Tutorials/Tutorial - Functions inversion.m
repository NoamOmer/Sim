clc
clear all
close all

% How to invert numerical functions in MATLAB
%
% Suppose you have a vector yy as a function of xx, meaning
% the command plot(xx,yy) is meaningful. For example, if the
% function f(x)=x^3 can be represented on the interval [-5 5] by:
xx=[-5:0.01:5];
yy=xx.^3;
plot(xx,yy);
title('y(x)')

% To invert it and obtain x(y), we will use the command
% interp1, syntax: YI = INTERP1(X,Y,XI). 
% Let xy denote the vector of x as a function of y, and let
% yyx denote the y-axis vector, the "new x-vector".
% We set yyx as follows: first find the min and max of y(x)

ymin = min(yy)
ymax = max(yy)

% Next set the resolution of our new vector yyx
N=100

% Now create it:
yyx=linspace(ymin,ymax,N);

% And now we interpolate xy along it:
xy=interp1(yy,xx,yyx);

% What we have actually done is said:
% The vector-pair yy-xx can also be viewed as xx-yy,
% in which xx is the y-values of a function x(y), given at
% NON-EQUIDISTANT points yy. We have simply interpolated
% x(y) at the EQUIDISTANT points yyx.

% Finally, just to check that it works:

figure
plot(yyx,xy)
title('x(y)');

% Viola!
