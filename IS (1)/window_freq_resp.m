% Input parameters
%   SW  : [kHz] : Spectral Width
%   SR  : [kHz] : Spectral Range
%   Tp  : [ms]  : Excitation time
%   ESP :       : Edge slope parameter

% Internal Parameters
%   N    :          : Number of points in the frequency response spectrum. Keeps dPhi/dw < PI.
%   df   : [kHz]    : Frequency axis resolution
%   Tp   : [ms]     : Excitation duration
%   freq : [kHz]    : 
%   R    : [kHz/ms] : Rate
function [windowed_freq_resp] = window_freq_resp(SW,SR,Tp,ESP,DEBUG_FLAG)

N    = max(5000,round(2*SW*(Tp*(SW/SR))));  %The term in the bracket is the extended Tp that keeps the same rate throughout the total SW
df   = SW/N;%1/20/Tp;
freq = -SW/2:df:SW/2-df;
R    = SR/Tp;
windowed_freq_resp = exp(-2*pi*1i*(freq).^2/2/R).*exp(-(1.95*abs(freq/SR)).^ESP);
windowed_freq_resp = windowed_freq_resp/max(abs(windowed_freq_resp));

if (DEBUG_FLAG >= 1)
	figure; hold;
	plot(freq, (abs(windowed_freq_resp)));
	title('Frequency Response');
	xlabel('\omega [kHz]');
end;

return;

