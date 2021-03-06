
function plot_post_acq_phase_1DUF(z_axis_,DNu0,OmegaE,te_of_z,Ga_of_ta,context);

set_context;

% Phi_e_total that was achieved
O_of_z      = -OmegaE;  % O(te(z)) = O(z)
dte_of_z    = ([(te_of_z(2) - te_of_z(1)) (te_of_z(2:end) - te_of_z(1:(end-1)))]);
Phi_RF_of_z = cumsum(O_of_z.*dte_of_z);               % [none]
Phi_e       = Phi_RF_of_z - OmegaE.*(Tp - te_of_z);   % [none]

% Phi_a
tmp1  = gammaHz*dta*cumsum(Ga_of_ta);
Phi_a = -tmp1(end)*z_axis_ - DNu0*Ta;

Phi_tot = Phi_e + Phi_a;
dPhi_tot_dz_numr = numeric_derivation(Phi_tot,z_axis_);

figure; hold;
plot(z_axis_,dPhi_tot_dz_numr*2*pi,'-');
title('d\phi total/dz @ t_a=T_a'); xlabel('Z-axis [cm]'); ylabel('d\Phi total/dz [rad]'); grid; set_gca;
legend({'d\phi/dz (achieved - num)'},'Location','Best');

return;

