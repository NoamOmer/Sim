% Input parameters
% ----------------
% RF_shape        : Shape of RF pulse          : [Rect_P, Sinc_P ...etc]
% RF_rot_phi      : Desired angle of rotation  : [rad]
% RF_phi_0        : Initial phase of RF pulse  : [rad]
% omega_cs        : Chemical shift frequency   : [rad]

% Internal parameters
% -------------------
% B_RF_amp_lab(t) : amplitude of RF pulse : lab frame : x-axis only --> 2 exponents --> 1 exponent --> cos+sin
% B_RF_amp_rot(t) : amplitude of RF pulse : rotating frame
% phi_RF_rot(t)   : phase of RF pulse     : rotating frame

% Output parameters
% -----------------
% B_RF_rot(t)     : 2D mat - [x(t),y(t),z(t)]
% B_RF_t(t)       : 1D vec - Time axis of RF pulse

function [B_RF_rot, B_RF_t] = gen_RF_pulse(pulse_prop, RF_rot_phi, RF_phi_0, dt, omega_cs)
set_globals;

B_RF_amp_lab = (2*pi)*B1MaxampLabHz/gamma_T;        % [T]       : omega[rad/sec] / gamma_T[rad/(sec*Tesla)]
B_RF_amp_rot = (1/2)*B_RF_amp_lab;                  % [T]       : half the amp of lab field
omega_RF     = gamma_T*B_RF_amp_rot;                % [rad/sec]

if (pulse_prop.type == Rect_P)
	[B_RF_t, pulse] = get_rectangular_pulse(RF_rot_phi, omega_RF, dt);
	B_RF_amp_rot = pulse * B_RF_amp_rot;
	phi_RF_rot   = (omega_0 + omega_cs)*B_RF_t + RF_phi_0 - omega_rot*B_RF_t;
elseif (pulse_prop.type == Sinc_P)
% 	[B_RF_t, pulse] = get_sinc_pulse(RF_rot_phi, omega_RF, dt, pulse_prop.n_lobes,pulse_prop.Tp);
	[B_RF_t, pulse] = get_sinc_pulse2(dt, pulse_prop.n_lobes,pulse_prop.Tp,pulse_prop.Hanning_Win_Flag);
	B_RF_amp_rot = pulse * B_RF_amp_rot;
	phi_RF_rot   = zeros(size(pulse));
% 	phi_RF_rot   = (omega_0 + omega_cs)*B_RF_t + RF_phi_0 - omega_rot*B_RF_t;
else
	error('gen_RF_pulse: Unknown RF shape (%d)', RF_shape);
end;

if (RH_flag)
	B_RF_rot(:,1) = transpose(B_RF_amp_rot .* cos(+phi_RF_rot));  % x
	B_RF_rot(:,2) = transpose(B_RF_amp_rot .* sin(+phi_RF_rot));  % y
else
	B_RF_rot(:,1) = transpose(B_RF_amp_rot .* cos(+phi_RF_rot));  % x
	B_RF_rot(:,2) = transpose(B_RF_amp_rot .* sin(-phi_RF_rot));  % y
end;
	B_RF_rot(:,3) = transpose(zeros(1,length(pulse)));            % z

if (pulse_prop.type == Sinc_P) && (pulse_prop.calib == 1)
	B_tmp = B_RF_rot(:,1) +1i*B_RF_rot(:,2);
	F = sum(B_tmp*gamma_T/(2*pi))*dt;
	flip_angle_deg  = RF_rot_phi*180/pi;
	flip_angle_cycl = (flip_angle_deg / 360);
	B_tmp = B_tmp * (flip_angle_cycl/F);
	F = sum(B_tmp*gamma_T/(2*pi))*dt;
	B_RF_rot(:,1) = real(B_tmp);
	B_RF_rot(:,2) = imag(B_tmp);
end;
	
if (DEBUG_FLAG >= 4)
	if (RH_flag)
		B1_RF_ = B_RF_amp_rot.*exp(+1i*phi_RF_rot);
	else
		B1_RF_ = B_RF_amp_rot.*exp(-1i*phi_RF_rot);
	end;
	figure;
	plot(B_RF_t*1E+3,1E-3*B_RF_amp_rot*gamma_T/2/pi,'b.-');
	title('RF(t) Amplitude');  xlabel('Time [ms]');  ylabel('B_1 Amp (kHz)');    set_gca;
	figure;
	plot(B_RF_t*1E+3,phi_RF_rot,'b.-');
	title('RF(t) Phase');      xlabel('Time [ms]');  ylabel('B_1 Phase [rad]');  set_gca;

	figure;
	plot(abs(fftshift(fft(B1_RF_))),'b.-');
	title('RF(\omega) Amplitude');  xlabel('none');  ylabel('B_1');  set_gca;
end;

return;

