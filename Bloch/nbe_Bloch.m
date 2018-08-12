
function [t, M_final, B_eff_rot] = nbe_Bloch(omega_cs, RF_shape, RF_rot_phi, RF_phi_0, plot_flag)

Bloch_set_globals;

B_lab = B0_lab + [0,0,(omega_cs/gamma_T)];
B_rot = B_lab  - [0, 0, (omega_rot/gamma_T)];

[B_RF_rot, RF_t] = gen_RF_pulse(RF_shape, RF_rot_phi, RF_phi_0, omega_cs, 'Bloch_set_globals');
B_RF_rot(end+1:length(t),1:3) = zeros(length(t) - length(B_RF_rot),3);               % pad RF pulse with zeroes
B_eff_rot = B_RF_rot + (ones(1,length(t))') * B_rot;    %'

n_pulses = 1;
for idx = 1:n_pulses
	M_final = apply_pulse_t(B_eff_rot, M0);
end;

if plot_flag
%	figure; plot(t,sqrt(B_rot(:,1).^2 + B_rot(:,2).^2),'.-'); title('Bxy Vs. t'); grid;

	t2 = 0:dt:n_pulses*(T+dt);
	figure;
	plot(t2,M_final(:,1),'r.-', t2,M_final(:,2),'b.-', t2,M_final(:,3),'k.-');
	legend({'Mx','My','Mz'});  grid;
end;

