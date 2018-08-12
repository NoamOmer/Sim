
function plot_pre_acq_phase_1DUF_RHR(z_axis_,Zia,Zfa,alpha0,alpha1,alpha2,DNu0,dDNu0_dz,INT_DNu0,OmegaE,te_of_z,Ge,O_of_te,Phi_RF_of_te,Tp,Ta);
set_globals;
global Ge;

% ==1==   Phi_e, according to which we design our RF.
%         Note that in some cases we havea PI pulse inside the sequence which reverses
%         the phase. We therefore design the RF to give a -phase which shows in the next
%         plots as opposite to the achieved / theoretical phase.
Phi_e     =   alpha2*(z_axis_.^2) + alpha1*z_axis_ + alpha0 - (Ta/(Zfa-Zia)) * (DNu0    .*(z_axis_ - Zia) - 2*INT_DNu0);
dPhi_e_dz = 2*alpha2*z_axis_      + alpha1                  - (Ta/(Zfa-Zia)) * (dDNu0_dz.*(z_axis_ - Zia) - DNu0      );
% dPhi_e_dz = numeric_derivation(Phi_e,z_axis_);

% ==2==   Phi_total that was achieved
O_of_z      = OmegaE;  % O(te(z)) = O(z)       % [Hz]
dte_of_z    = ([(te_of_z(2) - te_of_z(1)) (te_of_z(2:end) - te_of_z(1:(end-1)))]);
Phi_RF_of_z = cumsum(O_of_z.*dte_of_z);        % [none]

figure; hold; plot(1:length(z_axis_),O_of_z*1E-3,'-', 1:length(z_axis_),O_of_te*1E-3,'g--');
title('O(z) & O(t_e)'); xlabel('none'); ylabel('O [kHz]'); grid; set_gca;
legend({'O(z)','O(t_e)'},'Location','Best');

figure; hold; plot(1:length(z_axis_),Phi_RF_of_z*2*pi,'-', 1:length(z_axis_),Phi_RF_of_te,'g--');
title('\Phi_R_F(z) & \Phi_R_F(t_e)'); xlabel('none'); ylabel('\Phi_R_F [rad]'); grid; set_gca;
legend({'\Phi_R_F(z)','\Phi_R_F(t_e)'},'Location','Best');

Phi_tot = Phi_RF_of_z - (pi/2)/(2*pi) + OmegaE.*(Tp - te_of_z);
dPhi_tot_dz_numr = numeric_derivation(Phi_tot,z_axis_);
dPhi_tot_dz_ana  = (gammaHz*Ge + dDNu0_dz) .* (Tp - te_of_z);              % analytic calculation

% ==3==   Phi_total: theoretical phase value post Pi-pulse
Oi_ = O_of_z(1);
Of_ = O_of_z(end);
DO_ = abs(Of_ - Oi_);
Phi_tot_theoretical = - (((gammaHz*Ge)^2)*Tp/(2*DO_)).*(z_axis_.^2) + ...    % Theor. val of post Pi/2-pulse phase
                          (gammaHz*Ge*Tp*Of_ / DO_ ) .* z_axis_     + ...
                          Tp*((-Oi_)^2) / (2*DO_);

figure; hold;
plot(z_axis_,2*pi*Phi_e              ,'g-' );
plot(z_axis_,2*pi*Phi_tot            ,'--' );
plot(z_axis_,2*pi*Phi_tot_theoretical,'m-.');
title('\Phi total(z_e) @ t_a=0'); xlabel('z_e [cm]'); ylabel('\Phi total(z_e) [rad]'); grid; set_gca;
legend({'\phi_e desired'     ,...
        '\phi_e achieved'    ,...
        '\phi_e theoretical'},...
        'Location','Best');

figure; hold;
plot(z_axis_,2*pi*dPhi_e_dz  ,'g-');
plot(z_axis_,2*pi*dPhi_tot_dz_numr,'r-');
plot(z_axis_,2*pi*dPhi_tot_dz_ana ,'m--');
title('d\Phi total/dz @ t_a=0'); xlabel('Z-axis [cm]'); ylabel('d\Phi total/dz [rad/cm]'); grid; set_gca;
legend({'d\phi_e/dz (desired)'           ,...
        'd\phi_e/dz achieved (numeric)'  ,...
        'd\phi_e/dz achieved (analytic)'},...
        'Location','Best');

% dbstop in plot_pre_acq_phase_1DUF_RHR.m at 53
return;

