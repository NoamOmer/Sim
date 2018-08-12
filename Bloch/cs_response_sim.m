
clear all; close all; pack;
Bloch_set_globals;

omega_cs_arr = 2*pi*(-1E+3:0.01E+3:1E+3);   % [rad]
RF_shape   = Rect_P;
RF_rot_phi = pi/2;
plot_flag = 0;

M_final = [];
for idx = 1:length(omega_cs_arr)
	[t, M, B_eff_rot] = nbe_Bloch(omega_cs_arr(idx),RF_shape,RF_rot_phi,0,plot_flag);
	M_final(end+1,:) = M(end,:);
end;

freq_axis = 1E-3*omega_cs_arr/(2*pi);	% [kHz]
if (RF_shape == Rect_P)
	% Fig1: X,Y&Z
	figure; hold; grid;
	plot(freq_axis,M_final(:,1),'r- ');  plot(freq_axis,M_final(:,2),'b- ');  plot(freq_axis,M_final(:,3),'k.-');
	xlabel('Chemical shift [kHz]');    ylabel('Final magnetization ');	  legend({'M_x','M_y','M_z'});
	title('Final Magnetization (rotating frame)  Vs.  Chemical Shift');
	set_gca;
else
	% Fig1:Z  Fig2:X&Y
	figure; hold; grid;
	plot(freq_axis,M_final(:,3),'k.- '); title('Final M_z (rotating frame) Vs. Chemical Shift'); grid;
	set_gca;

% 	figure; hold; grid;
% 	plot(freq_axis,M_fin(:,1),'r- ');  plot(freq_axis,M_fin(:,2),'b- ');
% 	xlabel('Chemical shift [kHz]');    ylabel('Final magnetization ');    legend({'M_x','M_y'});
% 	title('Final M_x & M_y (rotating frame)  Vs.  Chemical Shift');
% 	set_gca;
end;

