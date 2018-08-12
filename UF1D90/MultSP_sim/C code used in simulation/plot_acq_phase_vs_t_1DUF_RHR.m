
% rows: z_axis;  columns: t
function plot_acq_phase_vs_t_1DUF_RHR(z_axis_,DNu0,dDNu0_dz,OmegaE,Tp,te_of_z,ta,dta_,Ga_of_ta,Ge);
set_globals;
global Ge;

% ----------------------------------------
% Phi(z)
% ----------------------------------------
% Phi_e_total that was achieved
O_of_z      = +OmegaE;                                                               % [Hz] O(te(z)) = O(z)
dte_of_z    = ([(te_of_z(2) - te_of_z(1)) (te_of_z(2:end) - te_of_z(1:(end-1)))]);   % [sec]
Phi_RF_of_z = cumsum(O_of_z.*dte_of_z);                                              % [none]
Phi_e       = Phi_RF_of_z - (pi/2)/(2*pi) + OmegaE.*(Tp - te_of_z);                  % [none]
Phi_e       = ones(length(ta),1)*Phi_e;                                              % [none]

% Phi_a
Phi_a = +gammaHz*transpose(cumsum(Ga_of_ta*dta_))*z_axis_ + transpose(ta)*DNu0;

% ------------------------------
% Total phase during acquisition
% ------------------------------
Phi_tot = Phi_e + Phi_a;    clear Phi_e Phi_a;

fh1 = figure; hold;
min_z = [];
min_phi = [];
for idx = [1:round(length(ta)/30):length(ta)]
	figure(fh1); hold on;
	plot(z_axis_,Phi_tot(idx,:),'k-','LineWidth',1);
	[min_phi(end+1),min_z(end+1)] = min(Phi_tot(idx,:));
end;
figure(fh1); hold on;
plot(z_axis_(min_z),min_phi,'rd-','LineWidth',1);
title('\phi(t,z)'); xlabel('Z-axis [cm]'); ylabel('Phase [none]'); grid; set_gca;
set(get(gca,'title' ),'FontWeight','BOLD','FontSize',24, 'FontName', 'Times New Roman');
set(get(gca,'ylabel'),'FontWeight','BOLD','FontSize',18, 'FontName', 'Times New Roman');
set(get(gca,'xlabel'),'FontWeight','BOLD','FontSize',18, 'FontName', 'Times New Roman');

clear O_of_z dte_of_z Phi_RF_of_z Phi_e Phi_a Phi_tot fh1 min_z min_phi idx; pack;

% ----------------------------------------
% dPhi/dz
% ----------------------------------------
% dPhi_e_dz that was achieved
dw_dz     = gammaHz*Ge + dDNu0_dz;                                                                    % [Hz/cm]
dPhi_e_dz = dw_dz .* (Tp - te_of_z);                                                                  % [1/cm]
dPhi_e_dz = ones(length(ta),1)*dPhi_e_dz;                                                             % [1/cm]

% dPhi_a_dz
dPhi_a_dz = gammaHz*transpose(cumsum(Ga_of_ta*dta_))*ones(1,length(z_axis_)) + transpose(ta)*dDNu0_dz; % [1/cm]

% Total phase during acquisition
dPhi_tot_dz = dPhi_e_dz + dPhi_a_dz;    clear dPhi_e_dz dPhi_a_dz;

fh2 = figure; hold;
for idx = [1:round(length(ta)/30):length(ta)]
	figure(fh2); hold on;
	plot(z_axis_,dPhi_tot_dz(idx,:),'k-','LineWidth',1);
end;
figure(fh2); hold on;
title('d\phi/dz(t,z)'); xlabel('Z-axis [cm]'); ylabel('Phase spatial derivative [none]'); grid; set_gca;

return;

