% Evovlves the input magnetization according the input magnetic field and acquire a signal through time.
% This function loops only on the temporal axis where in each temporal iteration we compute the effect
% of the current magnetic field on the entire sample as a bulk. I.e., no loop is preformed on the spatial axis.

% Note: This function assumes that Bx==0, By==0 & Bz~=0. I.e., the precession is solely in XY plane.
function [M_final, Sig_xy] = evolve_Mxy_n_acq(M_init_, B_eff_rot, z_grad, acquire_flag, context);

set_context;

M_final = M_init_;
Sig_xy = 0;

for t_idx = 1:length(B_eff_rot)   % Temporal forloop

	B_vec = B_eff_rot(t_idx,3)*ones(1,length(z_axis)) + z_grad*z_axis;   % Not very efficient... can be done outside the forloop
	rot_angle= gamma_T*B_vec*(t_idx*dt);
	
	new_M = (M_init_(:,1) + j*M_init_(:,2)) .* exp(-j*rot_angle');  %'

	M_final(:,1) = real(new_M);
	M_final(:,2) = imag(new_M);

	if (acquire_flag)
		M_xy = M_final(:,1) + j*M_final(:,2);  % M_xy = M_xy(z) at current t_idx
		Sig_xy(t_idx) = (1/sqrt(2*pi))*sum(M_xy * dz);
	end;
end;

return;

