
% 1D MRI Simulation

clear all; close all; pack;
MRI_1D_set_globals;

sample = exp(-1E+54*(((z_axis) * dz).^18));
% sample = 10 + 5*sin(4*pi*(1/Lz)*z_axis);

len_z = length(z_axis);
M_init = sample' * M0;   %'
B_lab = B0_lab + [0,0,(omega_cs/gamma_T)];
B_rot = B_lab  - [0, 0, (omega_rot/gamma_T)];

% ------------------
% Apply 90y RF pulse
% ------------------
[B_RF_rot, B_RF_t] = gen_RF_pulse(Rect_P, pi/2, 3*pi/2, 0, 0, 'MRI_1D_set_globals');
cur_grad = [0, 0, 0];
B_eff_rot = B_RF_rot + (ones(1,length(B_RF_t))') * B_rot;    %'
dB0z = zeros(1,length(z_axis));
[M, dummy] = evolve_M_n_acq(M_init, B_eff_rot, dB0z, cur_grad, dt, z_axis, 0, 0, 0);
% figure;
% plot(z_axis,M(:,1),'r.-', z_axis,M(:,2),'b>', z_axis,M(:,3),'k+-');
% title('M Vs. z - After 90y pulse');   legend({'M\_x','M\_y','M\_z'});

% ------------------------------------------------
% Apply negative gradient for duration of T_grad/2
% ------------------------------------------------
cur_grad = -GRAD(3);
B_grad = z_axis' * cur_grad;   %'
B_eff_rot = (ones(1,T_grad/(2*dt))') * B_rot;    %'
[M, dummy] = evolve_Mxy_n_acq(M, B_eff_rot, cur_grad, 0, 'MRI_1D_set_globals');

% Verification
% Mx_ver = M_init(:,3) .* cos(-gamma_T*cur_grad*z_axis*(T_grad/2))';   %'
% My_ver = M_init(:,3) .* sin(-gamma_T*cur_grad*z_axis*(T_grad/2))';   %'
% figure;
% plot(z_axis,M(:,1),'r.-', z_axis,Mx_ver,'m^');
% title('M_X Vs. z - After negative gradient for T\_grad/2');   legend({'M_X','M_X-Ver'});
% figure;
% plot(z_axis,M(:,2),'b.-', z_axis,My_ver,'c>');
% title('M_Y Vs. z - After negative gradient for T\_grad/2');   legend({'M_Y','M_Y-Ver'});

% -----------------------------------------------------------------
% Apply positive gradient for duration of T_grad and acquire signal
% -----------------------------------------------------------------
cur_grad = +GRAD(3);
B_grad = z_axis' * (GRAD);   %'
B_eff_rot = (ones(1,round(T_grad/dt))') * B_rot;    %'

[M, Sig_xy] = evolve_Mxy_n_acq(M, B_eff_rot, cur_grad, 1, 'MRI_1D_set_globals');

% Verification
% Mx_ver = M_init(:,3) .* cos(-gamma_T*GRAD(3)*z_axis*(T_grad/2))';   %'
% My_ver = M_init(:,3) .* sin(-gamma_T*GRAD(3)*z_axis*(T_grad/2))';   %'
% figure;
% plot(z_axis,M(:,1),'r.-', z_axis,Mx_ver,'m^-');
% title('M_X Vs. z - After positive gradient for T\_grad');   legend({'M_X','M_X-Ver'});
% figure;
% plot(z_axis,M(:,2),'b.-', z_axis,My_ver,'c>-');
% title('M_Y Vs. z - After positive gradient for T\_grad');   legend({'M_Y','M_Y-Ver'});

% Spatial reconstruction
% ----------------------
t_grad = 0:dt:(T_grad-dt);
k_t = gamma_T * GRAD(3)* ((t_grad - T_grad/2)');   %'
for idx = 1:length(z_axis)
	r_M(idx) = sum(Sig_xy .* exp(-j*k_t*z_axis(idx))' * dk);   %'
end;

figure; hold;
plot(z_axis,sample,'c.-'); plot(z_axis,abs(r_M),'m.-');
ylabel('Magnetization'); xlabel('z [m]');
title('Sample & Reconstructed sample');
legend({'Sample','Reconstructed sample'});

figure; hold; plot(abs(Sig_xy),'.-');
ylabel('FID'); xlabel('k [1/m]'); title('FID Vs k(t)');

