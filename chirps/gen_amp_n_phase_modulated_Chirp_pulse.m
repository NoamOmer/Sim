% Generate an Ampliture & Phase modulated RF Chirp pulse.

%  Parameters:
%  - - - - - -
%  te            :  [sec]      :  Time axis
%  R             :  [Hz/sec]   :  Rate dO/dte
%  Phi_RF        :  [rad]      :  RF phase
%  B1ampLabHz    :  [Hz]       :  RF field lab amplitude

%  Output:
%  - - - -
%  Complex (oscillating) exponential denoting the magnetic field in the rotating frame XY plane.
%  chirp_rot     :  [T]
%  phi_chirp_rot :  [rad]

function [chirp_rot, phi_chirp_rot] = gen_amp_n_phase_modulated_Chirp_pulse(te_,R_,Phi_RF_,context);
set_context;

% Time-Modulate RF amplitude - Since both Phi_RF and R are derived from O (integral and derivative respectively)
% we assume that the Phase Phi_RF already has the information concerning the sign of R.
% It is therefore ok to take the ABS of R and modulate the RF amplitude accordingly.
B1ampRotHz    = (1/2)*B1ampCoeff*sqrt(abs(R_));   % [Hz]
amp_chirp_rot = (2*pi*B1ampRotHz)/gamma_T;        % [T] = [rad/sec] / [rad/(sec*Tesla)]

% Time-Modulate RF phase
phi_chirp_rot = (omega_0 + omega_cs)*te_ + Phi_RF_ - omega_rot*te_;

if (DEBUG_FLAG >= 2)
	figure; hold;
    subplot(3,1,1); plot(te_*1E+3, R_*1E-6,'b-');
    title('Chirp Rate R(t)'); xlabel('Time [ms]'); ylabel('[kHz/ms]'); grid; set_gca;
    subplot(3,1,2); plot(te_*1E+3, amp_chirp_rot*gamma_T*1E-3/(2*pi),'r');
    title('Chirp Amplitude'); xlabel('Time [ms]'); ylabel('[kHz]');    grid; set_gca;
    subplot(3,1,3); plot(te_*1E+3, phi_chirp_rot,'m');
    title('Chirp Phase');     xlabel('Time [ms]'); ylabel('[rad]');    grid; set_gca;
end;

% RF pulse
if (RH_flag)
	chirp_rot = amp_chirp_rot .* exp(+j*phi_chirp_rot);
else
	chirp_rot = amp_chirp_rot .* exp(-j*phi_chirp_rot);
end;

return;

