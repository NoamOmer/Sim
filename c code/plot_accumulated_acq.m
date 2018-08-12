% Plot the signal from the entire sample, acquired at a specific point in time
% The signal is plotted Vs. the sample spatial axis,
% showing the contribution from each point along the smaple

% Mx0   :   vector   :   x-coordinate of the magnetization vector of the previous time epoch
% My0   :   vector   :   y-coordinate of the magnetization vector of the previous time epoch
% dt    :   double   :   Duration of previous time epoch
function plot_accumulated_acq(acc_real_vec,acc_imag_vec,Mx,My,Mz,Mx0,My0,dt,n_pixels,context)
set_context;
global pixel_idx;

if isempty(pixel_idx)
	pixel_idx = 0;
else
	pixel_idx = pixel_idx + 1;
end;
if (DEBUG_FLAG < 5)    return;    end;
Mxy  = Mx  + i*My;
Mxy0 = Mx0 + i*My0;
[dummy,low_idx]  = min(abs(sample_z_axis - z_axis(1)  ));
[dummy,high_idx] = min(abs(sample_z_axis - z_axis(end)));

% plot the magnetizaion vector
figure(104); hold on;
subplot(5,5,pixel_idx+1); hold on;
plot(sample_z_axis,abs(Mxy));
axis([-1.1 1.1 0 1.1]);
set(gca,'xticklabel',[],'yticklabel',[]);

% Calculate the phase and its derivatives Vs. z
Phi = phase(Mxy);                                      % [rad]
dPhi_dz   = numeric_derivation(Phi,sample_z_axis);     % [rad/cm]
d2Phi_dz2 = numeric_derivation(dPhi_dz,sample_z_axis); % [rad/cm^2]

% Calculate the current stationary point
[dPho_dz_val,dPho_dz_idx] = min(abs(dPhi_dz));
z_stationary = sample_z_axis(dPho_dz_idx);

% Calculate the current pixel size
pixel_size = abs(sqrt(2*pi/d2Phi_dz2(dPho_dz_idx)));   % [cm]

% figure(99);
% subplot(1,3,1);  hold on;  plot(z_stationary,dPho_dz_val,'.-');
% title('d\Phi/dz     Vs. Stationary point');   ylabel('d\Phi/dz [rad/cm]');   xlabel('z-axis [cm]');
% subplot(1,3,2);  hold on;  plot(z_stationary,pixel_size ,'.-');
% title('Pixel size   Vs. Stationary point');   ylabel('\Deltaz [cm]');        xlabel('z-axis [cm]');
% subplot(1,3,3);  hold on;  plot(1E+3*Ta*pixel_idx/n_pixels,z_stationary  ,'.-');
% title('Z-stationary Vs. t_a');                ylabel('z-stationary [cm]');   xlabel('Acq. Time [ms]');

% Shrink the vectors from sample to excitation size
dPhi_dz   = dPhi_dz(low_idx:high_idx);
d2Phi_dz2 = d2Phi_dz2(low_idx:high_idx);

% Calculate current temporal derivative of the phase
dPhi_dt = phase(Mxy0.*conj(Mxy))/(dt*2*pi);

acc_sig_vs_z = acc_real_vec(2:end) + i*acc_imag_vec(2:end);  % first cell in acc array is a dummy

% disp(sprintf('Pixel %4.0f(%3.0f). Size = %3.4d [cm]. z=%3.3f',pixel_idx,n_pixels,pixel_size,z_stationary));
% figure(100); hold on;
% subplot(3,2,1);  hold on;  plot(sample_z_axis, abs(Mxy));          xlabel('z-axis [cm]'); ylabel('M_{XY}(z)');
% title(['      M_{XY}(z)  Vs. z ' sprintf('(pix# %3.0f(%3.0f))',pixel_idx,n_pixels)]);
% 
% subplot(3,2,2);  hold on;  plot(sample_z_axis, dPhi_dt);           xlabel('z-axis [cm]'); ylabel('d\Phi/dt [rad/sec]');
% title(['    d\Phi/dt(z)  Vs. z ' sprintf('(pix# %3.0f(%3.0f))',pixel_idx,n_pixels)]);
% 
% subplot(3,2,3);  hold on;  plot(sample_z_axis, Phi);               xlabel('z-axis [cm]'); ylabel('\Phi_{M}(z) [rad]');
% title(['     \Phi_{M}(z) Vs. z ' sprintf('(pix# %3.0f(%3.0f))',pixel_idx,n_pixels)]);
% 
% subplot(3,2,4);  hold on;  plot(z_axis       , dPhi_dz);           xlabel('z-axis [cm]'); ylabel('d\Phi_{M}/dz [rad/cm]');
% title(['    d\Phi/dz (z) Vs. z ' sprintf('(pix# %3.0f(%3.0f))',pixel_idx,n_pixels)]);     axis([-1 1 -400 400]);
% 
% subplot(3,2,5);  hold on;  plot(sample_z_axis, abs(acc_sig_vs_z)); xlabel('z-axis [cm]'); ylabel('acc sig');
% title(['Accumulated Sig  Vs. z ' sprintf('(pix# %3.0f(%3.0f))',pixel_idx,n_pixels)]);     axis([-2 2 0 1020]);
% 
% subplot(3,2,6);  hold on;  plot(z_axis       , d2Phi_dz2);         xlabel('z-axis [cm]'); ylabel('d^2\Phi/dz^2 [rad/cm^2]');
% title(['d^2\Phi/dz^2 (z) Vs. z ' sprintf('(pix# %3.0f(%3.0f)) (Pix sz = %3.3f[cm])',pixel_idx,n_pixels,pixel_size)]);

% ------------------------------------------------------------------------------------------
% Calculate the current signal, given a uniformly excited sample of ones (S(x) = const. = 1)
% ------------------------------------------------------------------------------------------
% Use cumsum to see progression of signal accumulation
if 0 && (mod(pixel_idx+1,4) == 0)
	figure(101+pixel_idx); hold on;
	sig_uniform = cumsum(exp(i*Phi));
	sig         = cumsum(Mxy.*exp(i*Phi));
	zs = Zfe - Lz*(pixel_idx/n_pixels);
	str = sprintf('Zstat=%3.3f',zs);
	subplot(4,1,1); plot(sample_z_axis,   abs(sig_uniform)); xlabel('z-axis [cm]'); ylabel('Sig(z)'); title(['Uniform Sig(z) - Abs   ' str]);
	hold on; plot(zs,linspace(min(abs(sig_uniform))  ,max(abs(sig_uniform))  ,200),'r');
	subplot(4,1,2); plot(sample_z_axis, phase(sig_uniform)); xlabel('z-axis [cm]'); ylabel('Sig(z)'); title(['Uniform Sig(z) - Phase ' str]);
	hold on; plot(zs,linspace(min(phase(sig_uniform)),max(phase(sig_uniform)),200),'r');
	subplot(4,1,3); plot(sample_z_axis, abs(Mxy)          ); xlabel('z-axis [cm]'); ylabel('Sig(z)'); title(['Mxy(z) - Abs           ' str]);
	hold on; plot(zs,linspace(min(abs(Mxy))          ,max(abs(Mxy))          ,200),'r');
	subplot(4,1,4); plot(sample_z_axis, abs(sig)          ); xlabel('z-axis [cm]'); ylabel('Sig(z)'); title(['Sig(z) - Abs           ' str]);
	hold on; plot(zs,linspace(min(abs(sig))          ,max(abs(sig))          ,200),'r');
end;

% Use sum to see only final integral result
% global sum_idx;
% if (isempty(sum_idx))
% 	sum_idx = 1000;
% else
% 	sum_idx = sum_idx - 1;
% end;
sig_uniform = sum(exp(i*Phi));
sig         = sum(Mxy.*exp(i*Phi));
zs = Zfe - Lz*(pixel_idx/n_pixels);
str = sprintf('Zstat=%3.3f',zs);
figure(101);
subplot(3,1,1); hold on; plot(zs     ,   abs(sig_uniform),'.-'); xlabel('z-axis [cm]'); ylabel('Sig(z)'); title(['Uniform Sig(z) - Abs   ' str]);
subplot(3,1,2); hold on; plot(zs     , phase(sig_uniform),'.-'); xlabel('z-axis [cm]'); ylabel('Sig(z)'); title(['Uniform Sig(z) - Phase ' str]);
subplot(3,1,3); hold on; plot(zs     , abs(sig)          ,'.-'); xlabel('z-axis [cm]'); ylabel('Sig(z)'); title(['Sig(z) - Abs           ' str]);

% figure(102);
% subplot(3,1,1); hold on; plot(sum_idx,   abs(sig_uniform),'.-'); xlabel('z-axis [cm]'); ylabel('Sig(z)'); title(['Uniform Sig(z) - Abs   ' str]);
% subplot(3,1,2); hold on; plot(sum_idx, phase(sig_uniform),'.-'); xlabel('z-axis [cm]'); ylabel('Sig(z)'); title(['Uniform Sig(z) - Phase ' str]);
% subplot(3,1,3); hold on; plot(sum_idx, abs(sig)          ,'.-'); xlabel('z-axis [cm]'); ylabel('Sig(z)'); title(['Sig(z) - Abs           ' str]);

return;

