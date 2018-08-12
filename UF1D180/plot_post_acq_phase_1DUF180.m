
function plot_post_acq_phase_1DUF180(z_axis_,DNu0,OmegaE,Gpr,Ga_of_ta,context);

set_context;

% Phi_e_total that was achieved
O_of_z      = -OmegaE;  % O(te(z)) = O(z)
dte_of_z    = ([(te_of_z(2) - te_of_z(1)) (te_of_z(2:end) - te_of_z(1:(end-1)))]);
Phi_RF_of_z = cumsum(O_of_z.*dte_of_z);
Phi_e       = 2*Phi_RF_of_z - OmegaE.*(Tp - 2*te_of_z) - gammaHz*Gpr*Tpr*z_axis_;

% Phi_a
tmp1  = gammaHz*cumsum(Ga_of_ta*dta);
Phi_a = -tmp1(end)*z_axis_ - DNu0*Ta;

Phi_tot = Phi_e + Phi_a;
dPhi_tot_dz_numr = numeric_derivation(Phi_tot,z_axis_);

figure; hold;
plot(z_axis_,dPhi_tot_dz_numr ,'.-');
title('d\phi total/dz @ t_a=T_a'); xlabel('Z-axis [cm]'); ylabel('d\Phi total/dz'); grid; set_gca;
legend({'d\phi/dz (achieved - num)'},'Location','Best');

return;

