
function verify_inhom_fix_1D180(context);

set_context;
global Ge;

% ---------------------------------------------------------------------------------------------------------
% Verify the excitation range (this is an approximation since 'SR' should also include the B0 inhomo. term)
% ---------------------------------------------------------------------------------------------------------
SR = 1E-3*abs(gammaHz*Ge*Lz);                               % [1/sec]      Sample frequency range
if (abs(DOkHz) < SR)
	warning('verify_inhom_fix_1D180: Excitation range (%d[kHz]) is smaller than sample range (%d[kHz])',DOkHz,SR);
end;

% --------------------------------------------------
% Verify that the temporal resolution is high enough
% --------------------------------------------------
max_dte_required = abs(0.5*(1/(DOkHz*1E+3)));  % 0.5 factor, just for being on the safe side
if (dte > max_dte_required)
	error('verify_inhom_fix_1D180: Temporal resolution is too low (%d). Required resolution < %d',dte,max_dte_required);
end;

% ----------------------------------------------------------------------------------------
% Verify that the spatial resolution is small enough so as to reflect real-life experiment
% ----------------------------------------------------------------------------------------
if (exist('SIMULATION'))
	P_dDNu0_dz = polyder(P_DNu0);
	dDNu0_dz   = polyval(P_dDNu0_dz ,z_axis);

	if (RH_flag)
		dPhi_dz = +(gammaHz*Ge + dDNu0_dz).*(Tp - 2*te_of_z) + gammaHz*Gpr*Tpr;  % [1/cm]
	else
		dPhi_dz = -(gammaHz*Ge + dDNu0_dz).*(Tp - 2*te_of_z) - gammaHz*Gpr*Tpr;  % [1/cm]
	end;


	max_dz = (1/2)*min(abs(1./dPhi_dz));  % the factor 1/2 is just to be "on the safe side"
	if (dz > max_dz)
		error('verify_inhom_fix_1D180: Spatial resolution is too low (dz=%d, max_dz=%d',dz,max_dz);
	end;
end;

% ---------------------------------------------------
% Verify that the gradient strengths are not too high
% ---------------------------------------------------
if ((abs(Ge) > 50) || (abs(Gpr) > 50) || (max(abs(Ga_of_ta)) > 50))
	h = msgbox(sprintf('Warning: Ge=%3.1f, Gpr=%3.1f, max(abs(Ga))=%3.1f',Ge,Gpr,max(Ga_of_ta)));
	uiwait(h);
end;

