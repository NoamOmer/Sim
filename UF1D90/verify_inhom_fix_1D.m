
function verify_inhom_fix_1D(context)
set_context;
global Ge;

% ---------------------------------------------------------------------------------------------------------
% Verify the excitation range (this is an approximation since 'SR' should also include the B0 inhomo. term)
% ---------------------------------------------------------------------------------------------------------
SR  = 1E-3*abs(gammaHz*Ge(3)*Lz);              % [kHz]      Sample frequency range
if (abs(DOkHz) < SR)
	disp(sprintf('WARNING: Excitation range (%3.3f[kHz]) is smaller than sample range (%3.3f[kHz])',abs(DOkHz),SR));
end;

% --------------------------------------------------
% Verify that the temporal resolution is high enough
% --------------------------------------------------
%max_dte_required = abs(0.75*(1/(DOkHz*1E+3)));  % 0.75 factor, just for being on the safe side
max_dte_required = abs(0.90*(1/(DOkHz*1E+3)));   % 0.90 factor for experimental purposes
if (dte > max_dte_required)
	error('verify_inhom_fix_1D: Temporal resolution is too low (%3.1d). Required resolution < %3.1d',dte,max_dte_required);
end;

% ----------------------------------------------------------------------------------------
% Verify that the spatial resolution is small enough so as to reflect real-life experiment
% ----------------------------------------------------------------------------------------
if (exist('SIMULATION','var') && (SIMULATION == 1))
	P_dDNu0_dz = polyder(P_DNu0);
	dDNu0_dz   = polyval(P_dDNu0_dz ,z_axis);        clear P_dDNu0_dz;
	if (DesignPulse_Flag)
		if (RH_flag)
			dw_dz     = gammaHz*Ge(3) + dDNu0_dz;                            % [Hz/cm]
			dPhi_e_dz = dw_dz .* (Tp - te_of_z);         clear dw_dz;     % [1/cm]
			dPhi_e_dz = ones(length(ta),1)*dPhi_e_dz;                     % [1/cm]

			dPhi_a_dz = gammaHz*transpose(cumsum(Ga_of_ta*dta))*ones(1,length(z_axis)) + transpose(ta)*dDNu0_dz; %#ok<NODEF> % [1/cm]
		else
			dw_dz     =-gammaHz*Ge(3) - dDNu0_dz;                            % [Hz/cm]
			dPhi_e_dz = dw_dz .* (Tp - te_of_z);         clear dw_dz;     % [1/cm]
			dPhi_e_dz = ones(length(ta),1)*dPhi_e_dz;                     % [1/cm]

			dPhi_a_dz =-gammaHz*transpose(cumsum(Ga_of_ta*dta))*ones(1,length(z_axis)) - transpose(ta)*dDNu0_dz; %#ok<NODEF> % [1/cm]
		end;
	else
		dw_dz     = gammaHz*Ge(3) + dDNu0_dz;                                % [Hz/cm]
		dPhi_e_dz = dw_dz .* linspace(Tp,0,length(dw_dz));  clear dw_dz;  % [1/cm]
		dPhi_e_dz = ones(length(ta),1)*dPhi_e_dz;                         % [1/cm]

		Ga_of_ta = Ga(3)*ones(1,length(ta));
		dPhi_a_dz = gammaHz*transpose(cumsum(Ga_of_ta*dta))*ones(1,length(dPhi_e_dz)) + transpose(ta)*dDNu0_dz; % [1/cm]
	end;
	dPhi_tot_dz = dPhi_e_dz + dPhi_a_dz;                 clear dPhi_e_dz dPhi_a_dz;
	max_dz = min(min(abs(1./dPhi_tot_dz)));                               % [cm]
	max_dz = (1/2)*max_dz;           % add factor 1/2 just to be "on the safe side"
	if (dz > max_dz)
		error('verify_inhom_fix_1D: Spatial resolution is too low (dz=%5.3d, max_dz=%5.3d',dz,max_dz);
	end;
end;

% ---------------------------------------------------
% Verify that the gradient strengths are not too high
% ---------------------------------------------------
if (abs(Ge(3)) > 40) || (max(abs(Ga_of_ta)) > 40)
	h = msgbox(sprintf('WARNING: Ge = %3.1f, max(abs(Ga))=%3.1f',Ge,max(Ga_of_ta)));
	uiwait(h);
end;

return;

