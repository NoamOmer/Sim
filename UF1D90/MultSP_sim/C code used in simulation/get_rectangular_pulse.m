
% Generate a rectangular pulse.
%   rot_phi  : Required rotation angle of the magnetization
%   omega_RF : Precession frequency of the magnetization (depends on the strength of the magnetic field)
%   dt       : Temporal resolution

function [out_t_axis, out_pulse] = get_rectangular_pulse(rot_phi, omega_RF, dt);

% -------------------------------------------------------------------------------
% Disable this code -- in some cases we create more than one pulse per simulation
% -------------------------------------------------------------------------------
% global rect_pulse
% global pulse_t_axis;
% 
% if ~isempty(rect_pulse)
% 	out_pulse  = rect_pulse;
% 	out_t_axis = pulse_t_axis;
% 	return;
% end;

pulse_T      = rot_phi / omega_RF;       % [sec] = [rad] / [rad/sec]
pulse_len    = ceil(pulse_T/dt);         % number of pulse points
pulse_t_axis = 0:dt:((pulse_len-1)*dt);
rect_pulse   = ones(1,pulse_len);

out_pulse  = rect_pulse;
out_t_axis = pulse_t_axis;

return;

