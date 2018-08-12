% Generate an RF Chirp pulse.

%  Output:
%  - - - -
%  Complex (oscillating) exponential denoting the magnetic field in the rotating frame XY plane.

function [chirp_rot, phi_chirp_rot] = gen_Chirp_pulse(te_,Oi_,R_,B1ampCoeff_,context)
set_context;

B1ampRotHz    = B1ampCoeff_ * sqrt(abs(R_));      % [Hz]
amp_chirp_rot = (2*pi*B1ampRotHz)/gamma_T;        % [T] = [rad/sec] / [rad/(sec*Tesla)]

phi_chirp_rot = (omega_0 + omega_CS)*te_ + 2*pi*Oi_*te_ + (1/2)*2*pi*R_*(te_.^2) - omega_rot*te_;   % [rad]

if (RH_flag)
    chirp_rot = amp_chirp_rot * exp(+1i*phi_chirp_rot);
else
    chirp_rot = amp_chirp_rot * exp(-1i*phi_chirp_rot);
end;    

return;

% % OLD amplitude calculation code
% amp_chirp_rot = (1/2)*(2*pi*B1ampLabHz)/gamma_T;   % [T]

