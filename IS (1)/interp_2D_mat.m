
%       'nearest' - nearest neighbor interpolation
%       'linear'  - bilinear interpolation
%       'spline'  - spline interpolation
%       'cubic'   - bicubic interpolation as long as the data is
%                   uniformly spaced, otherwise the same as 'spline'

function [interpolated_m] = interp_2D_mat(m,F1,F2,interp_method)

[sx,sy] = size(m);
x = linspace(1,sx,sx*F1);
y = linspace(1,sy,sy*F2);

[xi,yi] = meshgrid(y,x);
if (exist('interp_method','var'))
	interpolated_m  = interp2(m,xi,yi,interp_method);
else
	interpolated_m  = interp2(m,xi,yi);
end;

return;

