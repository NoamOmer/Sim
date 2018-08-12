
% Input parameters
% ----------------
% M_init      :   []      :    2D matrix
%                              rows - spatial location along sample
%                              cols - x,y & z coordinates
% B_eff_rot   :   [T]      :   1D vector : [Bx By Bz]
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
function [M_final] = evolve_M(M_init_,B_eff_rot,dB0z_,grad,T,z_axis,inhomo_f);
set_globals; % Don't set context so as not to interfere with internal parameters

grad    = grad*1E-4;                           % Conversion to [T/cm]
B_grad  = [grad(1)*transpose(z_axis), grad(2)*transpose(z_axis), grad(3)*transpose(z_axis)];
M_final = M_init_;

% h1 = figure; hold;
% Phi_M_pre = angle(M_init_(:,1) + i*M_init_(:,2));
% plot(z_axis,Phi_M_pre,'b.--');
% title('\phi(z) Pre- & post- evolution phase');  xlabel('z-axis [cm]');  ylabel('Phase [rad]');
% 
% h2 = figure; hold;
% Phi_M_pre = phase(M_init_(:,1) + i*M_init_(:,2));
% plot(z_axis,Phi_M_pre,'b.--');
% title('\phi(z) Pre- & post- evolution phase');  xlabel('z-axis [cm]');  ylabel('Phase [rad]');

% h3 = figure;
% st_idx = 1;

for ax_idx = 1:length(M_final)   % Spatial forloop
% 	if (ax_idx == (round(length(M_final)/2) - 3))
% 		disp('stop')
% 	end;
	B_vec = B_eff_rot + B_grad(ax_idx,:);      % B_rot + spatial loc grad + ...
	if (inhomo_f)
		B_vec = B_vec + [0,0,dB0z_(ax_idx)];    % ... + B0 Inhomo.
	end;
	B_vec_sz = sqrt(sum(B_vec.^2));
	if (B_vec_sz == 0)
		R_mat = [1 0 0; 0 1 0; 0 0 1];         % No magnetic field --> no precession
	else;
		B_unit_vec = B_vec / B_vec_sz;
		phi_rot = gamma_T*B_vec_sz*T;          % [rad] = [rad/(sec*T)] * [T] * [sec]
		ph(ax_idx) = phi_rot;
        if (RH_flag)
			R_mat = rot_RHR(B_unit_vec, phi_rot);
        else
			R_mat = rot(B_unit_vec, phi_rot);
		end;
	end;
	M_final(ax_idx,:) = transpose(R_mat*transpose(M_final(ax_idx,:)));
    
	% T1 & T2 relaxation
 	if (Relax_Flag)
 		M_final(ax_idx,:) = T1_relaxation(M_final(ax_idx,:), MM_init_(3), T, T1);
 		M_final(ax_idx,:) = T2_relaxation(M_final(ax_idx,:), T, T2);
 	end;

% 	if (mod(ax_idx,30) == 0)
% 		figure(h3); hold on;
% 		Phi_M_post1 = phase(M_final(:,1) + i*M_final(:,2));
% 		plot(z_axis,Phi_M_post1,styles{st_idx});
% 		st_idx = st_idx+1;
% 	end;	
end;

% figure(h1);
% Phi_M_post = angle(M_final(:,1) + i*M_final(:,2));
% plot(z_axis,Phi_M_post,'r--');
% legend({'pre','post'});
% set_gca;
% 
% figure(h2);
% Phi_M_post = phase(M_final(:,1) + i*M_final(:,2));
% plot(z_axis,Phi_M_post,'r--');
% legend({'pre','post'});
% set_gca;

% figure;plot(z_axis,ph,'.-');title('Evolution phase change Vs z');
return;

