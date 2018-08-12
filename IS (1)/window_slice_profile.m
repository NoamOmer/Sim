
% n   [none]  length of vector
% dz  [cm]
% z0  [cm]
% fov [cm]
function [v] = window_slice_profile(n,dz,z0,fov,plot_f)

center_col = round(n * (1/2 + z0/fov) );
col_wdth   = round(n * dz/fov) - 1;
indices    = round(center_col - col_wdth/2) : round(center_col + col_wdth/2);
v = zeros(1,n);
v(indices) = 1;
smoothF = col_wdth;
if (smoothF <= 0)
	smoothF = 1;
end;
v = smooth(v,smoothF);

if (plot_f)
	figure;
	plot(linspace(-fov/2,fov/2,n),v,'.-');
	title('Slice profile');
	xlabel('[cm]');
	axis([-fov/2,fov/2,0,1.2]);
	grid;
end;
v = v / max(v);

return;

