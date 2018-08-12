
% rows: z_axis;  columns: t
function plot_acq_phase_vs_t_1DUF180_RHR(z_axis_,DNu0,dDNu0_dz,OmegaE,Tp,te_of_z,Gpr,Tpr,ta,dta,Ga_of_ta);
set_globals;
global Ge;

% Phi_e
O_of_z      = OmegaE;                                                                % [Hz] O(te(z)) = O(z)
dte_of_z    = ([(te_of_z(2) - te_of_z(1)) (te_of_z(2:end) - te_of_z(1:(end-1)))]);   % [sec]
Phi_RF_of_z = cumsum(O_of_z .* dte_of_z);                                            % [none]
Phi_e       = 2*Phi_RF_of_z + OmegaE.*(Tp - 2*te_of_z) + gammaHz*Gpr*Tpr*z_axis_;    % [none]
Phi_e       = ones(length(ta),1)*Phi_e;                                              % [none]

% Phi_a
Phi_a = gammaHz*transpose(cumsum(Ga_of_ta*dta))*z_axis_ + transpose(ta)*DNu0;

% ------------------------------
% Total phase during acquisition
% ------------------------------
Phi_tot = Phi_e + Phi_a;       clear Phi_e Phi_a;

fh1 = figure; hold;
if (Phi_tot(1,end) < Phi_tot(1,1))
	max_flag = 0;
else
	max_flag = 1;
end;	

stat_z = [];        % stationary point location
stat_phi = [];
for idx = [1:round(length(ta)/30):length(ta)]
	figure(fh1); hold on;
	plot(z_axis_,Phi_tot(idx,:),'k-','LineWidth',1);
	if (max_flag)
		[stat_phi(end+1),stat_z(end+1)] = max(Phi_tot(idx,:));
	else
		[stat_phi(end+1),stat_z(end+1)] = min(Phi_tot(idx,:));
	end;
end;

figure(fh1); hold on;
plot(z_axis_(stat_z),stat_phi,'rd-','LineWidth',1);
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
dPhi_e_dz = dw_dz .* (Tp - 2*te_of_z) + gammaHz*Gpr*Tpr;                                              % [1/cm]
dPhi_e_dz = ones(length(ta),1) * dPhi_e_dz;                                                           % [1/cm]

% dPhi_a_dz
dPhi_a_dz = gammaHz*transpose(cumsum(Ga_of_ta*dta))*ones(1,length(z_axis_)) + transpose(ta)*dDNu0_dz; % [1/cm]

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

