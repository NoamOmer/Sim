% Design the excitation pulse and acquisition gradient for fixing B0 inhomogeneity in 1DUF
% Uses Left-Hand Spin Rotation
%
%  Output parameters:
%  ------------------
%  te        :  1D vector            :  [sec]    :   Time axis
%  O_of_te   :  1D vector (te)       :  [Hz]     :   RF pulse time-dependant phase (dPhiRF_dte)
%  Phi_RF    :  1D vector (te)       :  [rad]    :   RF phase
%  Ge        :     scalar            :  [G/cm]   :   Excitation gradient
%  ta        :  1D vector (za_axis)  :  [sec]    :   Acquisition temporal axis
%  Ga_of_ta  :  1D vector (za_axis)  :  [G/cm]   :   Acquisition gradient

%  External parameters:
%  --------------------
%  Lz        :     scalar            :  [cm]     :   Sample length

%  Internal parameters:
%  --------------------
%  z_axis_   :  1D vector            :  [cm]     :   Sample spatial axis

function [te,dte,O_of_te,Ge,ta,Ga_of_ta,R_of_te,chirp_rot,phi_chirp_rot] = design_inhom_fix_1DUF(context);
set_context;
% Set the external 'z_axis' into the loacl 'z_axis_'
z_axis_ = z_axis;

% 'Z' axis for acquisition - These parems has to be defined here since we would like to be able to 
% use the current function without an external context
Zie_ = z_axis_(1);
Zfe_ = z_axis_(end);
Zia = Zfe_;
Zfa = Zie_;
disp(sprintf('Zie=%d, Zfe=%d, Zia=%d, Zfa=%d',Zie_,Zfe_,Zia,Zfa));

% Interpolate the internal za_axis so that the Ga dwell time in the ta axis 
% calculation will match the global dwell-time value
n_Ga_pts  = Ta/GRADdt;
hi_res_za = linspace(Zia,Zfa,n_Ga_pts);

if (DEBUG_FLAG >= 2)
	figure; hold; plot(1:length(z_axis_),z_axis_,'b.-');
	title('Excitation spatial axes'); ylabel('[cm]');
 	legend({'z_e axis'},'Location','Best'); set_gca; grid;
end;

% ----------------------------------------------------------------
% Calculate the needed derivations (and integrals) of Delta_Nu0(z)
% ----------------------------------------------------------------
P_INT_DNu0  = polyint(P_DNu0);
P_dDNu0_dz  = polyder(P_DNu0);

DNu0      = polyval(P_DNu0     ,z_axis_);   % [Hz]       Polynomial approximation for B0 inhomogeneity
INT_DNu0  = polyval(P_INT_DNu0 ,z_axis_);   % [Hz*cm]    spatial integral
dDNu0_dz  = polyval(P_dDNu0_dz ,z_axis_);   % [Hz/cm]    1st spatial derivative
if (DEBUG_FLAG >= 2
	figure; hold; plot(z_axis_,DNu0,'b-');
	title('Polynomial Approximation of Field Inhomogeneity \Delta\omega (z)');
	xlabel('z_e axis [cm]'); ylabel('\Delta\Omega [Hz]'); set_gca; grid;
end;

% -------------------
% Excitation gradient
% -------------------
alpha2 = 1/(2*(dZa^2));                                % [1/cm^2]
alpha1 = -(Ta/(Zfa-Zia))*(DNu0(end)) - 2*alpha2*Zfe_;  % [1/cm]    We assume that DNu0 is spanned along Ze and not Za
alpha0 = 0;                                            % [none]

% Calculate the excitation gradient: Ge
Ge = -(1/(gammaHz*Tp)) * ( (Ta/Lz)*DNu0(1) + (Tp+Ta)*dDNu0_dz(1) + 2*alpha2*z_axis_(1) + alpha1 );
disp(sprintf('Ge = %3.3f',Ge));

% --------------------
% Excitation time axis
% --------------------
% Calculate the effective frequency of spins as a fuction of 'z' during excitation: Omega(z)
OmegaE     = gammaHz*Ge*z_axis_ + DNu0 + omega_cs;
dOmegaE_dz = gammaHz*Ge         + dDNu0_dz;
dz_dOmegaE = 1./dOmegaE_dz;

OmegaE_min = min(OmegaE);
OmegaE_max = max(OmegaE);
OmegaE_axis  = linspace(OmegaE_min,OmegaE_max,length(z_axis_));
ze_of_OmegaE = interp1(OmegaE,z_axis_,OmegaE_axis);
dze_dOmegaE_app = diff(ze_of_OmegaE)/(OmegaE_axis(2)-OmegaE_axis(1));  % Approximated value
if (DEBUG_FLAG >= 2)
	figure; hold on;
	subplot(2,3,1); plot(z_axis_,OmegaE,'.-',z_axis_,gammaHz*Ge*z_axis_,'r.-');
	xlabel('z-axis [cm]');             ylabel('[Hz]');
	title('\omega_E(z)');              grid;    set_gca;
 	legend({'\omega_0+\Delta\omega_0','\omega_0'},'Location','Best');
 	
	subplot(2,3,2); plot(z_axis_,dOmegaE_dz,'.-');
	xlabel('z-axis [cm]');             ylabel('[Hz/cm]');
	title('d\omega_E/dz');             grid;    set_gca;
	
	subplot(2,3,3); plot(z_axis_,dz_dOmegaE,'.-');
	xlabel('z-axis [cm]');             ylabel('[cm/Hz]');
	title('1/(d\omega_E/dz)');         grid;    set_gca;

	subplot(2,3,4); plot(OmegaE_axis,ze_of_OmegaE,'.-',gammaHz*Ge*z_axis_,z_axis_,'r.-');
	xlabel('[Hz]');                    ylabel('[cm]');
	title('z(\omega_E) (app)');        grid;    set_gca;
 	legend({'z(\omega_0+\Delta\omega_0)','z(\omega_0)'},'Location','Best');
 	
	subplot(2,3,5); plot(ze_of_OmegaE(2:end),1./dze_dOmegaE_app,'.-');
	xlabel('z-axis [cm]');             ylabel('[Hz/cm]');
	title('d\omega_E/dz (app)');       grid;    set_gca;
	
	subplot(2,3,6); plot(ze_of_OmegaE(2:end),dze_dOmegaE_app,'.-');
	xlabel('z-axis [cm]');             ylabel('[cm/Hz]');
	title('1/(d\omega_E/dz) (app)');   grid;    set_gca;
end;

% Calculate the relation between time and location during excitation: te(z)
te_of_z = dz_dOmegaE .* (2*alpha2*z_axis_ + alpha1 - (Ta/Lz)*(dDNu0_dz.*(z_axis_ - z_axis_(end)) - DNu0)) + Tp;
if (DEBUG_FLAG >= 2)
	figure; hold; title('Excitation time of z-slices (t_e(z))');
	plot(z_axis_,te_of_z*1E+3,'-'); plot([z_axis_(1) z_axis_(end)], [te_of_z(1) te_of_z(end)]*1E+3, 'k--');
	xlabel('z-axis [cm]'); ylabel('t_e [ms]'); grid; set_gca;
end;
if ((sum(te_of_z < 0) > 1) || (min(te_of_z) ~= te_of_z(1)) || (max(te_of_z) ~= te_of_z(end)))
	error('te(z) has illegal values. Exiting');
end;

% Calculate a linear te axis
te  = linspace(min(te_of_z),max(te_of_z),length(z_axis_));
dte = te(2) - te(1);

% ----------------
% Excitation phase
% ----------------
ze_of_te = interp1(te_of_z,z_axis_,te);  % 'ze' at equaly spaces points along time axis 'te_axis'.
if (DEBUG_FLAG >= 2)
	figure; hold; title('z-Slice location along excitation time axis (z_e(t_e))');
	plot(te*1E+3,ze_of_te,'-'); plot([te(1) te(end)]*1E+3, [ze_of_te(1) ze_of_te(end)], 'k--');
	xlabel('t_e [ms]'); ylabel('z_e-axis [cm]'); grid; set_gca;
end;

% Calculate the RF phase temporal derivation O(te).
% Note that we first need to calculate DNu0 as a function of ze (as opposed to DNu0(z))
DNu0_of_ze_of_te = polyval(P_DNu0,ze_of_te);                      % [Hz]
O_of_te = -(gammaHz*Ge*ze_of_te + DNu0_of_ze_of_te + omega_cs);   % [Hz]  Chirp frequency (te)
Oi = O_of_te(1);                                                  % [Hz]  Initial excitation frequency
Phi_RF_of_te = cumsum(2*pi*O_of_te*dte);                          % [rad]
% Suggestion: maybe interpolate O(te) in order to get a more accurate result of the integral
if (DEBUG_FLAG > 1)
	figure; hold; title('O(t_e) = d\Phi_R_F(t_e)/dt');
	plot(te*1E+3, 1E-3*O_of_te             ,'-'  , ...
	     te*1E+3,-1E-3*gammaHz*Ge*ze_of_te ,'r--', ...
	     te*1E+3,-1E-3*DNu0_of_ze_of_te    ,'m-.');
	xlabel('t_e [ms]'); ylabel('O(t_e) [kHz]'); grid; set_gca;
    legend({'O(t_e)','-\gamma*G_e*z_e(t_e)','-\Delta\nu(z_e(t_e))'},'Location','Best');

    figure; hold; title('RF phase term due to \omega  vis-a-vis  the basic RF freqency');
	plot(te*1E+3, -cumsum(2*pi*gammaHz*Ge*ze_of_te*dte) ,'-'  , ...
	     te*1E+3, -cumsum(2*pi*DNu0_of_ze_of_te*dte)    ,'r--');
	xlabel('t_e [ms]'); ylabel('[rad]'); grid; set_gca;
    legend({'-2*\pi*\int(gammaHz*Ge*ze_of_te)dt_e','-2*\pi\int(DNu0_of_ze_of_te)dt_e'},'Location','Best');
end;

% ----------------
% Excitation power
% ----------------
R_of_te = numeric_derivation(O_of_te,te);    % [1/sec^2]    Chirp rate (dO/dte)
% !!!  Note  !!!
% the following API takes omega_cs into consideration. This means that it is accounted
% for twice instead of once, since it is also explicitly added to O_of_te.
[chirp_rot, phi_chirp_rot] = gen_amp_n_phase_modulated_Chirp_pulse(te,R_of_te,Phi_RF_of_te,context); % [[T],[rad]]
% if (DEBUG_FLAG)
% 	figure; hold; title('\Phi_R_F (t_e) [rad]'); xlabel('Time [ms]'); ylabel('[none]');
% 	plot(te*1E+3, phi_chirp_rot           ,'-' , ...
% 	     te*1E+3, unwrap(angle(chirp_rot)),'g--');
% 	legend({'\phi_R_F','unwrap(angle(RF))'},'Location','Best'); grid; set_gca;
% end;

% ---------------------
% Acquisition time axis
% ---------------------
% Create a time axis for acquisition
ta   = Ta * (hi_res_za - Zia) / (Zfa - Zia);
dta_ = ta(2) - ta(1);
if (DEBUG_FLAG >= 2)
	figure; hold; plot(ta*1E+3,hi_res_za,'-');
	title('Acquisition temporal axis t_a(z_a)'); xlabel('t_a [ms]'); ylabel('z_a-axis [cm]'); grid; set_gca;
end;

% -------------------- %
% Acquisition gradient %
% -------------------- %
% First calculate dDNu0 as a function of za.
% Since there's a linear (one-to-one) relation between za & ta it is sufficient to calculate d_DNu0/d_za
dDNu0_dza = polyval(P_dDNu0_dz,hi_res_za);
Ga_of_ta = -(1/gammaHz) * (2*alpha2*(Lz/Ta) + dDNu0_dza); -- do not use Lz in this eq. Replace with (Zfa-Zia) or something similar
if (DEBUG_FLAG >= 2)
	figure; hold; plot(ta*1E+3,Ga_of_ta,'-');
	title('Acquisition gradient G_a(t_a)'); xlabel('t_a-axis [ms]'); ylabel('G_a [G/cm]'); grid; set_gca;
end;

% ---------------------------------------------
% DEBUG: Total acquisition phase @ ta=0 & ta=Ta
% ---------------------------------------------
if (DEBUG_FLAG >= 3)
	plot_pre_acq_phase_1DUF(z_axis_,alpha0,alpha1,alpha2,DNu0,INT_DNu0,OmegaE,te_of_z,Ge,dDNu0_dz,Phi_RF_of_te,context);
	plot_post_acq_phase_1DUF(z_axis_,DNu0,OmegaE,te_of_z,Ga_of_ta,context);
	plot_acq_phase_vs_t_1DUF(z_axis_,DNu0,OmegaE,Tp,te_of_z,ta,dta_,Ga_of_ta);
end;

% -----------------------------------------------------------
% Ugly: necessary for external (mainly verification) purpuses
% -----------------------------------------------------------
% 'te' as a function of EXTERNAL 'z' axis -
% te_of_z_axis = interp1(z_axis_  ,te_of_z,z_axis);    % <-- uses EXTERNAL z_axis

% 'ta' as a function of EXTERNAL 'z' axis -
ta_of_z = interp1(hi_res_za, ta, z_axis);    % <-- uses EXTERNAL z_axis

return;

