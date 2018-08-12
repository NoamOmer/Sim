
function [hout] = plot_M_n_Phase(hin,sample_z_axis,M_init,M,M_title_str,Phase_title_str,Phase_leg)

global DEBUG_FLAG;
if (DEBUG_FLAG < 2)
	hout = 0;
	return;
end;
set_globals;

cur_fh = figure; subplot(211);
Mxy   = sqrt(M(:,1).^2 + M(:,2).^2);
Mxyz  = sqrt(M(:,1).^2 + M(:,2).^2 + M(:,3).^2);
plot(sample_z_axis*10,Mxy,'b.',sample_z_axis*10,M(:,3),'k--',sample_z_axis*10,Mxyz,'g-',sample_z_axis*10,M_init(:,3),'c-.');
title(M_title_str);  xlabel('z-axis [mm]');  ylabel('Magnetization'); % grid;
legend({'M_x_y','M_z','M_x_y_z','Init M_z'},'Location','Best');  set_gca;
a=axis; axis([a(1) a(2) -1.1 1.1]);
grid;

if (hin ~= 0)
    figure(hin); hold on;
    hout = hin;
    style = 'b-.';
else
    hout = figure(cur_fh);
	subplot(212);
    style = 'k.';
end;
Phi_M = unwrap(angle(M(:,1) + 1i*M(:,2)));
plot(sample_z_axis*10,Phi_M,style);
title(Phase_title_str);  xlabel('z-axis [mm]');  ylabel('Phase(z)'); grid;
legend(Phase_leg,'Location','Best');  set_gca;

return;
