
% Straighten linear offset, along the t2 (PE) k-space echo position
function sheared_kmat = shear(kmat,x_axis,shearing_factor,Phi0)
xmat = fftshift(fft(fftshift(kmat,2),[],2),2);

t2_axis = linspace(0,10,size(kmat,1));
[t2_axis_mat,x_axis_mat] = ndgrid(t2_axis,x_axis);
shear_mat = t2_axis_mat .* x_axis_mat .* shearing_factor;

% Add an arbitrary scaling factor
shear_mat = shear_mat ./ 25;

new_xmat = xmat .* exp(-i*shear_mat) .* exp(-i*Phi0);
sheared_kmat = ifftshift(ifft(ifftshift(new_xmat ,2),[],2),2);

return;

