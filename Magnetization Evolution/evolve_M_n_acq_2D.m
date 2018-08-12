% 2D (XY) Right-handed magnetization evolution in the presence of a 2D gradient

% Input parameters:
%  M_init      : []     :  2D cell of 1D vectors  :  Initial magnetization
%                          {Nx,Ny}                   Cell rows(1:end) - sample x coordinate x_initial --> x_final
%                                                    Cell cols(1:end) - sample y coordinate y_initial --> y_final
%                                                    1D vector - x,y,z magnetization values
%  B_eff_rot   : [T]    :  2D matrix [Nt,3]       :  Effective time dependent Rotating frame magnetic field
%                                                    [Bx(t) By(t) Bz(t)]
%                                                    rows - [Bx,By,Bz] at time 't'
%                                                    cols - [Bi(t)] (i=x,y,z)
%  dB0_XY      : [T]    :  2D matrix              :  B0 field inhomogeneity (in the z-direction)
%                          [Nx,Ny]                   rows(1:end) - sample x coordinate x_initial --> x_final
%                                                    cols(1:end) - sample y coordinate y_initial --> y_final
%  Grad_XY     : [G/cm] :  1D vector              :  X & Y gradients [Gx,Gy]
%  dt_         : [sec]  :  Scalar                 :  Time axis interval
%  x_axis      : [cm]   :  1D vector [Nx]         :  X Spatial axis
%  y_axis      : [cm]   :  1D vector [Ny]         :  Y Spatial axis
%  acq_flag    :        :                         :  Enable / disable signal acquisition
%  inhomo_flag :        :                         :  Enable / disable B0 inhomogeneity
%  plot_flag   :        :                         :  Enable / disable B0 inhomogeneity
%
% Output parameters:
%  M_final     : []     :  2D cell of 1D vectors  :  Final magnetization (see M_init for parameter format)
%  Sig_xy      : []     :  1D vector [Nt]         :  FID Signal Vs. time

function [M_final,Sig_xy] = evolve_M_n_acq_2D(M_init,B_eff_rot,dB0_XY,Grad_XY,dt_,x_axis,y_axis,...
                                                acq_flag,inhomo_flag,plot_flag);
set_globals;                         % Don't set context so as not to interfere with internal parameters

M_final   = M_init;

Grad_XY   = Grad_XY*1E-4;                                                 % Conversion [G/cm] --> [T/cm]
B_grad_XY = [Grad_XY(1)*transpose(x_axis) Grad_XY(2)*transpose(y_axis)];  % [Bgrad_X(x), Bgrad_Y(y)]

[Nx,Ny]   = size(M_final);
Nt        = size(B_eff_rot,1);

dx        = x_axis(2) - x_axis(1);
dy        = y_axis(2) - y_axis(1);
if (dx ~= dy)
	error(sprintf('evolve_M_n_acq_2D: dx(%2.3f) should be equal to dy(%2.3f)',dx,dy));
end;

Sig_xy    = 0;
if (plot_flag)  fh = figure;  end;

% In order to save the 'if' statement time inside the forloop - zero the inhomogeneity here if needed
if (~inhomo_flag)
	dB0_XY = dB0_XY*0;
end;

for t_idx = 1:Nt   % Temporal forloop
	if (mod(t_idx,20) == 0)  disp(sprintf('t_idx = %d (%d)',t_idx,length(B_eff_rot)));  end;

	for x_idx = 1:Nx              % X-axis spatial forloop
		for y_idx = 1:Ny          % Y-axis spatial forloop

			% Create an [Bx,By,Bz] vector including the gradient's induced field @ specific t
			% Remember that all gradients induce a magnetic field in the z-direction.
			B_vec = B_eff_rot(t_idx,:) + [0, 0, B_grad_XY(x_idx,1) + B_grad_XY(y_idx,2)];

			% Add the inhomogeneity contribution. Inhomogneity of B0 field, pointing in the z-direction.
			B_vec = B_vec + [0, 0, dB0_XY(x_idx,y_idx)];

			% Calculate the size of the B vector
			B_vec_sz = sqrt(sum(B_vec.^2));

			% Calculate the current rotation matrix
		    if (B_vec_sz == 0)
				R_mat = [1 0 0; 0 1 0; 0 0 1];     % No magnetic field --> no precession
			else
				B_unit_vec = B_vec / B_vec_sz;
				phi_rot = gamma_T*B_vec_sz*dt_;    % [rad] = [rad/(sec*T)] * [T] * [sec]
				R_mat = rot_RHR(B_unit_vec, phi_rot);
			end;
			
			% Rotate the Magnetization vector
			M_final{x_idx,y_idx} = transpose( R_mat * transpose(M_final{x_idx,y_idx}) );

			% T1 & T2 relaxation - commented out to save looping time
			% if (Relax_Flag)
			% 	M_final{x_idx,y_idx} = T1_relaxation(M_final{x_idx,y_idx}, M0(3), dt_, T1);
			% 	M_final{x_idx,y_idx} = T2_relaxation(M_final{x_idx,y_idx}       , dt_, T2);
			% end;
			
		end;
	end;

	if (acq_flag)
		M_final_vec  = [M_final{:}];
		Mx_final_vec = M_final_vec(1:3:length(M_final_vec));
		My_final_vec = M_final_vec(2:3:length(M_final_vec));
		M_xy = Mx_final_vec + j*My_final_vec;
		Sig_xy(t_idx) = sum(M_xy * dx * dy);
	end;

	% Comment out this part -- needed only for heavy debugging
	if (plot_flag & (mod(t_idx,50) == 0))
		M_final_vec  = [M_final{:}];
		Mx_final_vec = M_final_vec(1:3:length(M_final_vec));
		My_final_vec = M_final_vec(2:3:length(M_final_vec));
		Mx_final_mat = reshape(Mx_final_vec,Nx,Ny);
		My_final_mat = reshape(My_final_vec,Nx,Ny);
		Mxy = sqrt(Mx_final_mat.^2 + My_final_mat.^2);

		figure(fh); imagesc(x_axis,y_axis,Mxy);
		title(sprintf('M_x_y(x,y) at t=%d[ms]',t_idx*dt_*1E+3)); xlabel('x [cm]'); ylabel('y [cm]');
		pause(2);
	end;

end;

return;

