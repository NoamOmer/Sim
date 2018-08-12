% Calculate the largest step in the signal phase.
% Assumptions:
%   1.  y_axis is the excitation axis.
%   2.  Acquisition is done in the reverse direction (from y_axis end to start)

% y_axis    [cm]
% alpha0    [none]
% alpha1    [1/(cm)]
% alpha2    [1/(cm^2)]
% Ga        [G/cm]
% Ta        [sec]
% Nacq      [none]
% SE_flag   [none]

function [status, dPhi_max] = SPEN_acq_phase_check_aliasing(y_axis,alpha0,alpha1,alpha2,Ga,Ta,Nacq,SE_flag)
gammaHz = 4257.4;    % [Hz/G]

if (Nacq ~= length(y_axis))
	disp('spatial / time axes mismatch');
	status = 0;
	return;
end;

phi_e = alpha2*(y_axis.^2) + alpha1*y_axis + alpha0;
if (SE_flag), phi_e = -phi_e;  end;

ta = linspace(0,Ta,Nacq);
k  = gammaHz*Ga*ta;
dt = Ta / Nacq;

  dPhi_max = 2*pi * ( (phi_e(end-1) + gammaHz*Ga*dt   *y_axis(end-1)) - (phi_e(end)                               ) );   % same but even simpler
% dPhi_max = 2*pi * ( (phi_e(end-1) + gammaHz*Ga*dt   *y_axis(end-1)) - (phi_e(end) + gammaHz*Ga*0    *y_axis(end)) );   % same but simpler
% dPhi_max = 2*pi * ( (phi_e(end-1) + gammaHz*Ga*ta(2)*y_axis(end-1)) - (phi_e(end) + gammaHz*Ga*ta(1)*y_axis(end)) );   % original

if (abs(dPhi_max) < pi)
	status = 1; % ok
else
	status = 0; % nok
end;

% figure; hold on;
% for idx = 1:1:Nacq
% plot(y_axis,2*pi*(phi_e + k(idx)*y_axis),'k.-');
% end;
% plot(y_axis,2*pi*(phi_e),'b.-');
% plot(y_axis,2*pi*(phi_e + k(end)*y_axis),'r.-');

return;




% 
% function [status] = SPEN_acq_phase_check_aliasing(y_axis,alpha0,alpha1,alpha2,Ga,Ta,Nacq,SE_flag)
% 
% if (Nacq ~= length(y_axis))
% 	disp('spatial / time axes mismatch');
% 	status = 0;
% 	return;
% end;
% 
% L = abs(y_axis(end) - y_axis(1));
% ya_start = y_axis(end);
% ya_end   = y_axis(1);
% 
% t_axis = (y_axis + ya_start) * Ta / L;
% 
% phi_e = alpha2*(ye_axis.^2) + alpha1*ye_axis + alpha0;
% if (SE_flag)
% 	phi_e = -phi_e;
% end;
% phi_a_atStationaryPoint = phi_e + gammaHz*Ga*t_axis.*ya_axis;
% 
% dPhi_max = phi_a_atStationaryPoint(2) - phi_a_atStationaryPoint(1);
% 
% function [status] = SPEN_acq_phase_check_aliasing(ye_axis,alpha0,alpha1,alpha2,Ga,Ta,Nacq,SE_flag)
% 
% if (Nacq ~= length(ye_axis))
% 	disp('spatial / time axes mismatch');
% 	status = 0;
% 	return;
% end;
% 
% ya_axis = ye_axis(end:-1:1);
% 
% L = abs(ye_axis(end) - ye_axis(1));
% ye_start = y_axis(1);
% ye_end   = y_axis(end);
% 
% t_axis = (y_axis + ye_start) * Ta / L;
% 
% phi_e = alpha2*(ye_axis.^2) + alpha1*ye_axis + alpha0;
% if (SE_flag)
% 	phi_e = -phi_e;
% end;
% phi_a_atStationaryPoint = phi_e + gammaHz*Ga*t_axis.*ya_axis;
% 
% dPhi_max = phi_a_atStationaryPoint(2) - phi_a_atStationaryPoint(1);
% 
