
% Generate a SINC - sin(x)/x - pulse.

function [out_t_axis, out_pulse] = get_sinc_pulse2(dt, n_lobes, Tp, Hanning_Win_Flag)

RF_T = Tp;

pulse_len = ceil(RF_T/dt);
pulse_t_axis = 0:dt:((pulse_len-1)*dt);

% dbstop in get_sinc_pulse2.m at 12
% sinc_pulse = sinc(linspace(-n_lobes/2,+n_lobes/2,pulse_len));
sinc_pulse = sin(pi*linspace(-n_lobes/2,+n_lobes/2,pulse_len)) ./ (pi*linspace(-n_lobes/2,+n_lobes/2,pulse_len));
if (Hanning_Win_Flag)
sinc_pulse = sinc_pulse .* transpose(hanning(length(sinc_pulse)));
end;

out_pulse  = sinc_pulse;
out_t_axis = pulse_t_axis;

return;

