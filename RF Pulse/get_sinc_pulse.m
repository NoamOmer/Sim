
% Generate a SINC - sin(x)/x - pulse.

function [out_t_axis, out_pulse] = get_sinc_pulse(phi_rot, omega_RF, dt, n_lobes, Tp)

% global sinc_pulse;
% global pulse_t_axis;
% 
% if ~isempty(sinc_pulse)
% 	out_pulse  = sinc_pulse;
% 	out_t_axis = pulse_t_axis;
% 	return;
% end;

RF_T = Tp; %phi_rot / omega_RF;

pulse_len = ceil(RF_T/dt);
pulse_t_axis = 0:dt:((pulse_len-1)*dt);

sinc_pulse = sinc(linspace(-n_lobes/2,+n_lobes/2,pulse_len));

% sinc_end = round(pulse_len/2);
% ratio = (pulse_len/15)/2;
% sinc_axis = (-sinc_end:1:sinc_end)/ratio;
% sinc_axis = sinc_axis(1:pulse_len);
% sinc_pulse = sinc(sinc_axis);

% out_pulse  = sinc_pulse*28;    % Factor '28' was set experimentally in order to match the power of Rect pulse
out_pulse  = sinc_pulse;
out_t_axis = pulse_t_axis;

return;

