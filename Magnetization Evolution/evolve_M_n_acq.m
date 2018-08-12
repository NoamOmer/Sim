 
% Input parameters
% ----------------
% M_init      :   []      :    2D matrix
%                              rows - spatial location along sample
%                              cols - x,y & z coordinates
% B_eff_rot   :   [T]      :   2D matrix : [Bx(t) By(t) Bz(t)]
%                              rows - temporal location along time axis
%                              cols - x,y & z coordinates
% grad        :   [G/cm]   :   1D vector : [Gx Gy Gz]
%                              Magnetic field gradient [T/m]. Assumption: constant gradient, G == [0 0 Gz(x,y,z)]
% dB0z        :   [T]      :   1D vector : dB0z(z)
%                              Magnetic field inhomogeneity along the spatial axis
%
% Output parameters
% -----------------
% M_final   : Magnetization after time evolving. Same structure as M_init.
% Sig_xy    : 1D vector - MRI Signal Vs. time
function [M_final,Sig_xy] = evolve_M_n_acq(M_init_,B_eff_rot,dB0z_,grad,dt_,z_axis,acquire_flag,inhomo_f,plot_f)
set_globals; % Don't set context so as not to collide with internal parameters

grad    = grad*1E-4;                                             % Conversion to [T/cm]
B_grad  = [grad(1)*transpose(z_axis), grad(2)*transpose(z_axis), grad(3)*transpose(z_axis)];
M_final = M_init_;
Sig_xy  = 0;
dz      = z_axis(2) - z_axis(1);
fh      = 0;
len_M     = length(M_final);

% In order to save the 'if' statement time inside the forloop - zero the inhomogeneity here if needed
if (~inhomo_f)
	dB0z_ = dB0z_*0;
end;

%WB=createWaitBar('init',0);
%createWaitBar('title',sprintf('Magnetization evolution (acq=%d)',acquire_flag),WB);
for t_idx = 1:length(B_eff_rot)   % Temporal forloop
%	createWaitBar('percent',t_idx/length(B_eff_rot)*100,WB);
%	if WB.error
%		error('evolve_M_n_acq1Dgrad_t: Error in Waitbar. Exiting.');
%	end;

	if (mod(t_idx,100) == 0)
		disp(sprintf('t_idx = %d (%d)',t_idx,length(B_eff_rot)));
	end;
	for ax_idx = 1:len_M          % Spatial forloop
		B_vec = B_eff_rot(t_idx,:) + B_grad(ax_idx,:);   % Current B_rot + spatial loc grad + ...
		B_vec = B_vec + [0,0,dB0z_(ax_idx)];             % ... + B0 Inhomo.
		B_vec_sz = sqrt(sum(B_vec.^2));
	    if (B_vec_sz == 0)
			R_mat = [1 0 0; 0 1 0; 0 0 1];               % No magnetic field --> no precession
		else
			B_unit_vec = B_vec / B_vec_sz;
			phi_rot = gamma_T*B_vec_sz*dt_;              % [rad] = [rad/(sec*T)] * [T] * [sec]
            if (RH_flag)
				R_mat = rot_RHR(B_unit_vec, phi_rot);
            else
				R_mat = rot(B_unit_vec, phi_rot);
			end;
		end;
		M_final(ax_idx,:) = transpose(R_mat*transpose(M_final(ax_idx,:)));
	    
		% T1 & T2 relaxation
		if (Relax_Flag)
			M_final(ax_idx,:) = T1_relaxation(M_final(ax_idx,:), M_init_(3), dt_, T1);
			M_final(ax_idx,:) = T2_relaxation(M_final(ax_idx,:), dt_, T2);
		end;
	end;

% Comment out this part -- needed only for serious debugging
%	if (plot_f & (mod(t_idx,100) == 0))
%		if (~fh) fh = figure; end;
%		Mxy = sqrt(M_final(:,1).^2 + M_final(:,2).^2);
%		figure(fh); plot(z_axis,Mxy,'b.-'); axis([1.2*min(z_axis) 1.2*max(z_axis) -2.1 2.1]);
%		title(sprintf('M_x_y(z) at t=%d[ms]',t_idx*dt_*1E+3)); xlabel('z [cm]'); ylabel('Magnetization');
%		pause(0.1E+0);
%	end;
	
	if (acquire_flag)
		M_xy = M_final(:,1) + j*M_final(:,2);  % M_xy = M_xy(z) at current t_idx
		Sig_xy(t_idx) = sum(M_xy * dz);
	end;
end;
%delete(WB.figure);

return;

