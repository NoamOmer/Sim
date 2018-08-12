% Design a selective soft pulse.

% Input parameters
% ----------------
% RF_shape        : Shape of RF pulse          : [Rect_P, Sinc_P ...etc]
% RF_Amp          : Amplitude of RF pulse      : [Hz]
% dt_             : Pulse temporal resolution  : [sec]
% center_freq     : Center frequency           : [Hz]
% freq_wdth       : Width of selective pulse   : [Hz]
% phi0            : Initial phase of RF pulse  : [rad]

% Internal parameters
% -------------------
% B_RF_amp_lab(t) : amplitude of RF pulse : lab frame : x-axis only --> 2 exponents --> 1 exponent --> cos+sin
% B_RF_amp_rot(t) : amplitude of RF pulse : rotating frame
% phi_RF_rot(t)   : phase of RF pulse     : rotating frame

% Output parameters
% -----------------
% B_RF_rot(t)     : 2D mat - [x(t),y(t),z(t)]
% B_RF_t(t)       : 1D vec - Time axis of RF pulse

! TBD !
function [B_RF_rot,B_RF_t] = gen_presat_pulse(RF_shape,flip_a,dt_,center_freq,freq_wdth,phi0,context);
set_context;

pulse_T    = 1/freq_wdth;                           % [sec]     : The selective bandwidth = 1 / pulse duration
omega_RF   = flip_a / pulse_T;                      % [rad/sec]

%B_RF_amp_lab = (2*pi)*RF_Amp/gamma_T;               % [T]       : omega[rad/sec] / gamma_T[rad/(sec*Tesla)]
%B_RF_amp_rot = (1/2) *B_RF_amp_lab;                 % [T]       : half the amp of lab field
%omega_RF = gamma_T*B_RF_amp_rot;                    % [rad/sec]


if (RF_shape == Rect_P)
	[B_RF_t, pulse] = get_rectangular_pulse(flip_angle, omega_RF, dt_);
elseif (RF_shape == Sinc_P)
	[B_RF_t, pulse] = get_sinc_pulse       (flip_angle, omega_RF, dt_);
else
	error('gen_RF_pulse: Unknown RF shape (%d)', RF_shape);
end;

B_RF_amp_rot = pulse * B_RF_amp_rot;
phi_RF_rot   = (omega_0 + omega_cs - omega_rot)*B_RF_t + center_freq*2*pi*B_RF_t + phi0;

if (RH_flag)
	B_RF_rot(:,1) = transpose(B_RF_amp_rot .* cos(+phi_RF_rot));
	B_RF_rot(:,2) = transpose(B_RF_amp_rot .* sin(+phi_RF_rot));
else
	B_RF_rot(:,1) = transpose(B_RF_amp_rot .* cos(+phi_RF_rot));
	B_RF_rot(:,2) = transpose(B_RF_amp_rot .* sin(-phi_RF_rot));
end;
	B_RF_rot(:,3) = transpose(zeros(1,length(pulse)));

if (DEBUG_FLAG)
	if (RH_flag)
		B1_RF_ = B_RF_amp_rot.*exp(+i*phi_RF_rot);
	else
		B1_RF_ = B_RF_amp_rot.*exp(-i*phi_RF_rot);
	end;
	figure;
	plot(B_RF_t*1E+3,1E-3*B_RF_amp_rot*gamma_T/2/pi,'b.-');
	title('Pre-Sat(t) Amplitude');  xlabel('Time [ms]');  ylabel('B_1 Amp (kHz)');    set_gca;
	figure;
	plot(B_RF_t*1E+3,phi_RF_rot,'b.-');
	title('Pre-Sat(t) Phase');      xlabel('Time [ms]');  ylabel('B_1 Phase [rad]');  set_gca;

	figure;
    df   = 1/B_RF_t(end);
    fmax = 1/(B_RF_t(2)-B_RF_t(1));
    f    = -fmax/2 : df : (fmax/2);
	plot(f,abs(fftshift(fft(B1_RF_))),'b.-');
	title('Pre-Sat(\omega) Amplitude');  xlabel('Frequency [Hz]');  ylabel('B_1');  set_gca;
end;

return;

