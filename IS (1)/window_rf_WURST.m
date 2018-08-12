% Input
%   rf   : vector : [???]  : Time dependant RF shape
%   t_rf : vector : [sec]  : Time axis
%   ESP  : scalar : [none] : Edge slope parameter
function [windowed_rf_amp] = window_rf_WURST(rf_amp,t_rf,ESP,DEBUG_FLAG,title_str)

Tp = t_rf(end);
RFwin = 1 - abs((sin(pi*(t_rf-Tp/2)/Tp)).^ESP);
windowed_rf_amp = rf_amp .* RFwin;

if (DEBUG_FLAG >= 1)
	figure; hold;
	plot(t_rf, (abs(rf_amp         )*1E+4)*(4.2574e+3*1E-3),'b.-');
	plot(t_rf, (abs(windowed_rf_amp)*1E+4)*(4.2574e+3*1E-3) ,'g.-');
	title('RF Pulse(t)');
	xlabel('t_{RF} [sec]');
	ylabel('RF Amp');
	legend({'pre-win','post-win'},'Location','Best');
end;

return;

