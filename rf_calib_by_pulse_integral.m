
% Calibrate RF according to the pulse integral and the required flip-angle
% b1     : complex : RF shape              [a.u.]
% angle  : real    : flip-angle            [deg]
% out_b1 : complex : calibrated RF shape   [T]
function out_b1 = rf_calib_by_pulse_integral(b1,dt,flip_angle_deg)
set_globals;

% Normalize input
b1_tmp = b1 / max(abs(b1));

% calculate normalized flip-angle for this pulse via:
% alpha = integral((B1[T] * gamma_T[rad/(sec*T] / [rad]) * dt[sec])
norm_b1_flip_angle  = sum(b1_tmp*gamma_T/(2*pi))*dt;

% convert flip angle to cycles
flip_angle_cycl = (flip_angle_deg / 360);

% calibrate b1
out_b1 = b1_tmp * (flip_angle_cycl/norm_b1_flip_angle);

% self-consistency check
% new_b1_flip_angle  = sum(out_b1*gamma_T/(2*pi))*dt;
% fprintf('required   flip-angle = %3.1f [deg]\n',flip_angle_deg       );
% fprintf('calibrated flip-angle = %3.1f [deg]\n',new_b1_flip_angle*360);

return;
