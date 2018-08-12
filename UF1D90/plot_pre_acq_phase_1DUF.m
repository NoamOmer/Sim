
function plot_pre_acq_phase_1DUF(z_axis_,alpha0,alpha1,alpha2,DNu0,INT_DNu0,OmegaE,te_of_z,Ge,dDNu0_dz,Phi_RF_of_te,context);

set_context;
global Ge;

% ==1==   Phi_e, according to which we design our RF & purge parameters
Phi_e = alpha2*(z_axis_.^2) + alpha1*z_axis_ + alpha0 - (Ta/Lz) * (DNu0.*(z_axis_ - z_axis_(end)) - 2*INT_DNu0);
dPhi_e_dz = numeric_derivation(Phi_e,z_axis_);

% ==2==   Phi_total that was achieved
O_of_z      = -OmegaE;  % O(te(z)) = O(z)       % [Hz]
dte_of_z    = ([(te_of_z(2) - te_of_z(1)) (te_of_z(2:end) - te_of_z(1:(end-1)))]);
Phi_RF_of_z = cumsum(O_of_z.*dte_of_z);        % [none]

figure; hold; plot(1:length(z_axis_),O_of_z*1E-3,'-', 1:length(z_axis_),O_of_te*1E-3,'g--');
title('O(z) & O(t_e)'); xlabel('none'); ylabel('O [kHz]'); grid; set_gca;
legend({'O(z)','O(t_e)'},'Location','Best');

figure; hold; plot(1:length(z_axis_),Phi_RF_of_z,'-', 1:length(z_axis_),Phi_RF_of_te/(2*pi),'g--');
title('\Phi_R_F(z) & \Phi_R_F(t_e)'); xlabel('none'); ylabel('\Phi_R_F [rad]'); grid; set_gca;
legend({'\Phi_R_F(z)','\Phi_R_F(t_e)'},'Location','Best');

Phi_tot = Phi_RF_of_z - OmegaE.*(Tp - te_of_z) - 0;     % '0': zero terms, from acq phase
dPhi_tot_dz_numr = numeric_derivation(Phi_tot,z_axis_);
dPhi_tot_dz_ana  = -(gammaHz*Ge + dDNu0_dz) .* (Tp - te_of_z);  % analytic calculation

% ==3==   Phi_total: theoretical phase value post Pi-pulse
Oi = O_of_z(1);
Of = O_of_z(end);
DO_Hz = abs(Of - Oi);
Phi_tot_theoretical = (((gammaHz*Ge)^2)*Tp/(2*DO_Hz)).*(z_axis_.^2) + ...    % Theor. val of post Pi/2-pulse phase
                      -(gammaHz*Ge*Tp*Of / DO_Hz ) .* z_axis_ - ...
                       Tp*((-Oi)^2) / (2*DO_Hz);

figure; hold;
plot(z_axis_,Phi_e  ,'g-');
plot(z_axis_,Phi_tot,'--');
plot(z_axis_,Phi_tot_theoretical,'m-.');
title('\Phi total(z_e) @ t_a=0'); xlabel('z_e [cm]'); ylabel('\Phi total(z_e)'); grid; set_gca;
legend({'\phi_e desired','\phi_e achieved','\phi_e theoretical'},'Location','Best');

figure; hold;
plot(z_axis_,dPhi_e_dz  ,'g-');
plot(z_axis_,dPhi_tot_dz_numr,'r-');
plot(z_axis_,dPhi_tot_dz_ana ,'m--');
title('d\Phi total/dz @ t_a=0'); xlabel('Z-axis [cm]'); ylabel('d\Phi total/dz'); grid; set_gca;
legend({'d\phi_e/dz (desired)','d\phi_e/dz (achieved - num)','d\phi_e/dz (achieved - ana)'},'Location','Best');

return;

