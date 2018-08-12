
function plot_pre_acq_phase_1DUF180(z_axis_,alpha0,alpha1,alpha2,DNu0,INT_DNu0,OmegaE,Gpr,Ge,dDNu0_dz,Ga_of_ta,Phi_RF_of_te,context);

set_context;

% == Phi_e ==  according to which we design our RF & purge parameters
Phi_e = alpha2*(z_axis_.^2) + alpha1*z_axis_ + alpha0 - (Ta/Lz) * (DNu0.*(z_axis_ - z_axis_(end)) - 2*INT_DNu0); % [none]
dPhi_e_dz1 = numeric_derivation(Phi_e,z_axis_); % [1/cm]
dPhi_e_dz2 = 2*alpha2*z_axis_ + alpha1 - (Ta/Lz) * (dDNu0_dz.*(z_axis_ - z_axis_(end)) - DNu0); % [1/cm]

% == Phi_RF ==
% according to desired Phi_e
Phi_RF_of_z1 = 0.5 * ( Phi_e + OmegaE.*(Tp - 2*te_of_z) + gammaHz*Gpr*Tpr*z_axis_ );

% according to integral of Omega(z)
O_of_z       = -OmegaE;  % [Hz] O(te(z)) = O(z)
dte_of_z     = ([(te_of_z(2) - te_of_z(1)) (te_of_z(2:end) - te_of_z(1:(end-1)))]);  % [sec]
% Unsuccessful - Try to add the constant value of Phi_RF_of_z1 to Phi_RF_of_z2
% P = polyfit(z_axis,Phi_RF_of_z1,100);
% P(end)
% Phi_RF_of_z2 = cumsum(O_of_z .* dte_of_z) + P(end);
Phi_RF_of_z2 = cumsum(O_of_z .* dte_of_z);                                           % [none]

figure; hold;
plot(te*1E+3,O_of_z,'-', te*1E+3,O_of_te,'g-',te*1E+3,numeric_derivation(Phi_RF_of_te/(2*pi),te),'k-');
title('O(z), O(t_e)   &   d\Phi_R_F(t_e)/dt_e (numeric)'); xlabel('[none]   &   Time [ms]'); ylabel('O [Hz]');
legend({'O(z)','O(t_e)','d\Phi_R_F/dt_e (num.)'},'Location','Best'); grid; set_gca;

figure; hold;
plot(1:length(z_axis_), Phi_RF_of_z2       ,'b-' ,...
     1:length(z_axis_), Phi_RF_of_te/(2*pi),'k-',...
     1:length(z_axis_), Phi_RF_of_z1       ,'r-');
title('\Phi_R_F(z) & \Phi_R_F(t_e)'); xlabel('none'); ylabel('\Phi_R_F [none]'); grid; set_gca;
legend({'\Phi_R_F(z) - Int(\omega(z))','\Phi_R_F(t_e) - Int(\omega(t_e))','\Phi_R_F(z) - analytical'},'Location','Best');

% == Phi_total ==  that was achieved - at ta=0
Phi_tot = 2*Phi_RF_of_z2 - OmegaE.*(Tp - 2*te_of_z) - gammaHz*Gpr*Tpr*z_axis_;     % [none]
dPhi_tot_dz_numr = numeric_derivation(Phi_tot,z_axis_);                            % [1/cm] numeric  calculation
dPhi_tot_dz_ana  = -(gammaHz*Ge + dDNu0_dz) .* (Tp - 2*te_of_z) - gammaHz*Gpr*Tpr; % [1/cm] analytic calculation

% == Phi_total ==  theoretical phase value post Pi-pulse
Oi = O_of_z(1);    % [Hz]
Of = O_of_z(end);  % [Hz]
DO_Hz= Of - Oi;    % [Hz]
Phi_tot_theoretical = -(((gammaHz*Ge)^2)*Tp/DO_Hz).*(z_axis_.^2) + ...          % Theor. val of post Pi-pulse phase
                        (gammaHz*Ge*Tp*(Of + Oi) / DO_Hz ) .* z_axis_ - ...     % Linear part should be = 0
                        Tp*(Oi^2) / DO_Hz;                                      % [none]

figure; hold;
plot(z_axis_,Phi_e  ,'b-');
plot(z_axis_,Phi_tot,'r-');
%	plot(z_axis_,Phi_tot_theoretical,'k-');
title('\Phi total(z_e) @ t_a=0'); xlabel('z_e [cm]'); ylabel('\Phi total(z_e)'); grid; set_gca;
%	legend({'\phi_e desired','\phi_e achieved','\phi_e theoretical'},'Location','Best');
legend({'\phi_e desired','\phi_e achieved'},'Location','Best');

figure; hold;
plot(z_axis_,Phi_e-Phi_tot  ,'-');
title('\Phi total(z_e) @ t_a=0  :  desired-achieved'); xlabel('z_e [cm]'); ylabel('\Delta\Phi total(z_e)'); grid; set_gca;

figure; hold;
plot(z_axis_,dPhi_e_dz1  ,'b-');
plot(z_axis_,dPhi_tot_dz_numr,'r-');
plot(z_axis_,dPhi_tot_dz_ana ,'k-');
title('d\Phi total/dz @ t_a=0'); xlabel('Z-axis [cm]'); ylabel('d\Phi total/dz'); grid; set_gca;
legend({'d\phi_e/dz (desired)','d\phi_e/dz (achieved - num)','d\phi_e/dz (achieved - ana)'},'Location','Best');

return;

