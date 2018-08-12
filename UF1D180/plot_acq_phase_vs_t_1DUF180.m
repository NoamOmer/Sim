% rows: z_axis;  columns: t

function plot_acq_phase_vs_t_1DUF180(z_axis_,DNu0,OmegaE,Tp,te_of_z,Gpr,Tpr,ta,dta,Ga_of_ta);
set_globals;
global F1;
global F2;
% zDecimationFactor = 5;
% tDecimationFactor = 100;
% z_axis_  = z_axis_ (1:zDecimationFactor:end);
% OmegaE   = OmegaE  (1:zDecimationFactor:end);
% te_of_z  = te_of_z (1:zDecimationFactor:end);
% DNu0     = DNu0    (1:zDecimationFactor:end);
% Ga_of_ta = Ga_of_ta(1:tDecimationFactor:end);
% ta       = ta      (1:tDecimationFactor:end);
% dta      = ta(2) - ta(1); % dta changes due to ta decimation

% -----------------------------
% Phi_e total that was achieved
% -----------------------------
O_of_z      = -OmegaE;  % O(te(z)) = O(z)
dte_of_z    = ([(te_of_z(2) - te_of_z(1)) (te_of_z(2:end) - te_of_z(1:(end-1)))]);
Phi_RF_of_z = cumsum(O_of_z.*dte_of_z);
Phi_e       = 2*Phi_RF_of_z - OmegaE.*(Tp - 2*te_of_z) - gammaHz*Gpr*Tpr*z_axis_;
Phi_e       = ones(length(ta),1)*Phi_e;

% -----
% Phi_a
% -----
Phi_a = -gammaHz*transpose(cumsum(Ga_of_ta*dta))*z_axis_ - transpose(ta)*DNu0;

% ------------------------------
% Total phase during acquisition
% ------------------------------
Phi_tot = Phi_e + Phi_a;

h1 = figure; hold;
h2 = figure; hold;
min_z = [];
min_phi = [];
F1_(1) = getframe;
F2_(1) = getframe;
for idx = [1:round(length(ta)/20):length(ta)]
	cur_phi_tot = Phi_tot(idx,:);
	cur_dphi_tot_dz = numeric_derivation(cur_phi_tot,z_axis_);
	figure(h1); hold on;
	plot(z_axis_,cur_phi_tot,'k-','LineWidth',1);
%		ax=axis;
%		axis([z_axis_(1), z_axis_(end),ax(3),ax(4)]);
	[min_phi(end+1),min_z(end+1)] = min(cur_phi_tot);
	plot(z_axis_(min_z),min_phi,'rv-','LineWidth',2);
    F1_(end+1)=getframe;

	figure(h2); hold on;
	plot(z_axis_,cur_dphi_tot_dz,'k-','LineWidth',1);
end;

h3 = figure;
F2_(1) = getframe;
for idx = [1:round(length(ta)/150):length(ta)]
	cur_phi_tot = Phi_tot(idx,:);
	figure(h3);
    ph = mod(cur_phi_tot*360/2/pi,360);
	plot(z_axis_,ph,'k-','LineWidth',1);
    F2_(end+1)=getframe;
end;
F1 = F1_;
F2 = F2_;
%     movie2avi(F1,'AcqPhase_vs_t','fps',20,'compression','none');
%     movie2avi(F1,'AcqPhase_vs_t2');
figure(h1);
title('\phi(t,z)'); xlabel('Z-axis [cm]'); ylabel('Phase'); grid; set_gca;
set(get(gca,'title' ),'FontWeight','BOLD','FontSize',24, 'FontName', 'Times New Roman');
set(get(gca,'ylabel'),'FontWeight','BOLD','FontSize',18, 'FontName', 'Times New Roman');
set(get(gca,'xlabel'),'FontWeight','BOLD','FontSize',18, 'FontName', 'Times New Roman');

figure(h2);
title('d\phi(t,z)/dz'); xlabel('Z-axis [cm]'); ylabel('Phase Derivative'); grid; set_gca;
set(get(gca,'title' ),'FontWeight','BOLD','FontSize',24, 'FontName', 'Times New Roman');
set(get(gca,'ylabel'),'FontWeight','BOLD','FontSize',18, 'FontName', 'Times New Roman');
set(get(gca,'xlabel'),'FontWeight','BOLD','FontSize',18, 'FontName', 'Times New Roman');

return;

