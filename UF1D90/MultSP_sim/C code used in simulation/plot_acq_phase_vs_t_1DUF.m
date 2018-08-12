
% rows: z_axis;  columns: t
function plot_acq_phase_vs_t_1DUF(z_axis_,DNu0,OmegaE,Tp,te_of_z,ta,dta_,Ga_of_ta);
set_globals;

% Phi_e_total that was achieved
O_of_z      = -OmegaE;  % O(te(z)) = O(z)
dte_of_z    = ([(te_of_z(2) - te_of_z(1)) (te_of_z(2:end) - te_of_z(1:(end-1)))]);
Phi_RF_of_z = cumsum(O_of_z.*dte_of_z);               % [none]
Phi_e       = Phi_RF_of_z - OmegaE.*(Tp - te_of_z);   % [none]
Phi_e       = ones(length(ta),1)*Phi_e;

% Phi_a
Phi_a = -gammaHz*dta_*transpose(cumsum(Ga_of_ta))*z_axis_ - transpose(ta)*DNu0;

% Total phase during acquisition
Phi_tot = Phi_e + Phi_a;

figure; hold;
min_z = [];
min_phi = [];
for idx = [1:round(length(ta)/20):length(ta)]
	plot(z_axis_,Phi_tot(idx,:),'k-','LineWidth',1);
	[min_phi(end+1),min_z(end+1)] = min(Phi_tot(idx,:));
end;
plot(z_axis_(min_z),min_phi,'rd-','LineWidth',1);
title('\phi(t,z)'); xlabel('Z-axis [cm]'); ylabel('Phase [none]'); grid; set_gca;
set(get(gca,'title' ),'FontWeight','BOLD','FontSize',24, 'FontName', 'Times New Roman');
set(get(gca,'ylabel'),'FontWeight','BOLD','FontSize',18, 'FontName', 'Times New Roman');
set(get(gca,'xlabel'),'FontWeight','BOLD','FontSize',18, 'FontName', 'Times New Roman');

return;


