
function plot_B0z_inhomo(B0z_lab, dB0z, z_axis);

set_globals;

figure; hold;
plot(100*dB0z./B0z_lab, z_axis, '.-');
title('B_0 Inhomogeneity - dB_0 (% of B_0)')
xlabel('dB_0 [% of B_0]'); ylabel('z-axis [cm]');
set_gca;

return;

