
%       'nearest' - nearest neighbor interpolation
%       'linear'  - bilinear interpolation
%       'spline'  - spline interpolation
%       'cubic'   - bicubic interpolation as long as the data is
%                   uniformly spaced, otherwise the same as 'spline'

function [resized_m] = resize_mat(m,sz1,sz2,interp_method)

[sx,sy] = size(m);
x = linspace(1,sx,sz1);
y = linspace(1,sy,sz2);

[xi,yi] = meshgrid(y,x);
if (exist('interp_method','var'))
	resized_m = interp2(m,xi,yi,interp_method);
else
	resized_m = interp2(m,xi,yi);
end;

return;
