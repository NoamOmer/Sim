% Design the excitation pulse and acquisition gradient for fixing B0 inhomogeneity in 2D Hybrid UF/EPI
% Ver11 -- Evolve magnetization according to the Right-Hand-Rule
% Ver12 -- Set Ge instead of dZa
function [te,dte,O_of_te,Ge,Gpr,Ga_of_ta,R_of_te,chirp_rot,phi_chirp_rot] = design_inhom_fix_UF_1D180_RHR(context);
set_context;

z_axis_ = z_axis;   % Set the external 'z_axis' into the loacl 'z_axis_'

% ------------------
% Acquisition z-axis
% ------------------
% 'Z' axis for acquisition - These parems has to be defined here since we would like to be able to 
% use the current function without an external context
Zie_ = z_axis_(1);
Zfe_ = z_axis_(end);
Zia  = Zfe_;
Zfa  = Zie_;
disp(sprintf('Zie=%3.3f, Zfe=%3.3f, Zia=%3.3f, Zfa=%3.3f',Zie_,Zfe_,Zia,Zfa));

% Interpolate the internal za_axis so that the Ga dwell time in the ta axis 
% calculation will match the global dwell-time value
n_Ga_pts  = Ta/GRADdt;
hi_res_za = linspace(Zia,Zfa,n_Ga_pts);                       % curent
% hi_res_za = [Zia:(dZa/Za_AXIS_INTERP_FACTOR):Zfa];          % old
% hi_res_za = z_axis(end):(z_axis(1)-z_axis(2)):z_axis(1);    % older

if (DEBUG_FLAG >= 2)
	figure; hold; plot(1:length(z_axis_),z_axis_,'b.-');
	title('Excitation spatial axis'); ylabel('z_e axis [cm]');
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
if (DEBUG_FLAG >= 2)
	figure; hold; plot(z_axis_,DNu0,'b.-');
	title('Polynomial Approximation of Field Inhomogeneity \Delta\omega(z)');
	xlabel('z_e axis [cm]'); ylabel('\Delta\Omega [Hz]'); set_gca; grid;
end;

% ----------------
% Purge parameters
% ----------------
Gpr    = (Tp / (gammaHz*Tpr)) * (gammaHz*Ge + dDNu0_dz(end));
alpha2 = (1/(2*(Zie_-Zia))) * ( (Ta/(Zfa-Zia)) * ( dDNu0_dz(1)*(Zie_-Zia) - DNu0(1) + DNu0(end) ) + ...
                                Tp * (gammaHz*Ge + dDNu0_dz(1)) + ...
                                gammaHz*Gpr*Tpr );                       % [1/cm^2]
alpha1 = -2*alpha2*Zia - (Ta/(Zfa-Zia))*DNu0(end);                       % [1/cm]
alpha0 = 0;                                                              % [none]
dZa    = sign(alpha2) * 1/sqrt(2*abs(alpha2));

% ------------------------------
% Validation
% ------------------------------
disp(sprintf('\n Tpr=%d\n Gpr=%d\n Tp=%d\n Ge=%d\n Ta=%d\n dZa=%d',Tpr,Gpr,Tp,Ge,Ta,dZa));
if (Map_Mult_Factor == 0)
	if (abs(Ge) - abs(alpha2*Lz/(gammaHz*Tp)) > abs(Ge/10))
 		error('Ge is not equal to alpha2*Lz/(gammaHz*Tp) for zero inhomogeneity');
	end;
end;

% --------------------
% Excitation time axis
% --------------------
% Calculate the effective frequency of spins as a fuction of 'z' during excitation: Omega(z)
OmegaE     = gammaHz*Ge*z_axis_ + DNu0 + omega_cs;
dOmegaE_dz = gammaHz*Ge         + dDNu0_dz;
dz_dOmegaE = 1./dOmegaE_dz;

OmegaE_axis     = linspace(min(OmegaE),max(OmegaE),length(z_axis_));
ze_of_OmegaE    = interp1(OmegaE,z_axis_,OmegaE_axis);
dze_dOmegaE_app = numeric_derivation(ze_of_OmegaE,OmegaE_axis); %diff(ze_of_OmegaE)/(OmegaE_axis(2)-OmegaE_axis(1));  % Approximated value
if (DEBUG_FLAG >= 2)
	figure; hold on;
	subplot(2,3,1); plot(z_axis_,1E-3*OmegaE,'.-',z_axis_,1E-3*gammaHz*Ge*z_axis_,'r.-');
	xlabel('z-axis [cm]');            ylabel('[kHz]');
	title('\omega_E(z)');             grid;    set_gca;
 	legend({'\omega_0+\Delta\omega_0','\omega_0'},'Location','Best');
 	
	subplot(2,3,2); plot(z_axis_,1E-3*dOmegaE_dz,'.-');
	xlabel('z-axis [cm]');            ylabel('[kHz/cm]');
	title('d\omega_E/dz');            grid;    set_gca;

	subplot(2,3,3); plot(z_axis_,1E+3*dz_dOmegaE,'.-');
	xlabel('z-axis [cm]');            ylabel('[cm/kHz]');
	title('1/(d\omega_E/dz)');        grid;    set_gca;

	subplot(2,3,4); plot(1E-3*OmegaE_axis,ze_of_OmegaE,'.-',1E-3*gammaHz*Ge*z_axis_,z_axis_,'r.-');
	xlabel('[kHz]'); ylabel('[cm]');
 	legend({'z(\omega_0+\Delta\omega_0)','z(\omega_0)'},'Location','Best');
 	title('z(\omega_E) (app)');       grid;    set_gca;

	subplot(2,3,5); plot(ze_of_OmegaE,1E-3*1./dze_dOmegaE_app,'.-');
	xlabel('z-axis [cm]');            ylabel('[kHz/cm]');
	title('d\omega_E/dz (app)');      grid;    set_gca;

	subplot(2,3,6); plot(ze_of_OmegaE,1E+3*dze_dOmegaE_app,'.-'); 
	xlabel('z-axis [cm]');            ylabel('[cm/kHz]');
	title('1/(d\omega_E/dz) (app)');  grid;    set_gca;
end;
% Calculate the relation between time and location during excitation: te(z)
te_of_z = 0.5 * Tp - ...
          0.5 * (1 ./ (gammaHz*Ge + dDNu0_dz)) .* ( 2*alpha2*z_axis_ + alpha1 - ...
                                                    (Ta/(Zfa-Zia)) * (dDNu0_dz.*(z_axis_ - Zia) - DNu0) - ...
                                                    gammaHz*Gpr*Tpr );
if (DEBUG_FLAG >= 2)
	figure; hold; title('Excitation time of z-slices (t_e(z))');
	plot(z_axis_,te_of_z*1E+3,'.-'); plot([z_axis_(1) z_axis_(end)], [te_of_z(1) te_of_z(end)]*1E+3, 'k--');
	xlabel('z-axis [cm]'); ylabel('t_e [ms]'); grid; set_gca;
end;
if ((sum(te_of_z < 0) > 2) || (min(te_of_z) ~= te_of_z(1)) || (max(te_of_z) ~= te_of_z(end)))
	figure; hold; title('Excitation time of z-slices (t_e(z))');
	plot(z_axis_,te_of_z*1E+3,'.-'); plot([z_axis_(1) z_axis_(end)], [te_of_z(1) te_of_z(end)]*1E+3, 'k--');
	xlabel('z-axis [cm]'); ylabel('t_e [ms]'); grid; set_gca;
	h = msgbox('te(z) has illegal values. Would you like to continue?');
	uiwait(h);
end;

% Calculate a linear te axis
te  = linspace(min(te_of_z),max(te_of_z),length(z_axis_));
dte = te(2) - te(1);

% ----------------
% Excitation phase
% ----------------
% Calculate ze(te) during excitation by reversing te(z)
ze_of_te = interp1(te_of_z,z_axis_,te);  % 'ze' at equaly spaced points along time axis 'te'.
if (DEBUG_FLAG >= 2)
	figure; hold; title('z-Slice location along excitation time axis (z_e(t_e))');
	plot(te*1E+3,ze_of_te,'-',te*1E+3,z_axis_,'g--');
	xlabel('t_e [ms]'); ylabel('z_e-axis [cm]'); grid; set_gca;
	legend({'z_e(t_e)','z_e axis'},'Location','Best');
end;

% Calculate O(te) - the temporal derivative of the RF phase
% Note that we first need to calculate DNu0(ze(te)) (as opposed to DNu0(z))
DNu0_of_ze_of_te = polyval(P_DNu0,ze_of_te);
O_of_te = (gammaHz*Ge*ze_of_te + DNu0_of_ze_of_te + omega_cs);         % [Hz]
Phi_RF_of_te = cumsum(2*pi*O_of_te*dte);                               % [rad]
% Suggestion: maybe interpolate O(te) in order to get a more accurate result of the integral
if (DEBUG_FLAG >= 2)
	figure; hold; title('O(t_e) = d\Phi_R_F(t_e)/dt');
	plot(te*1E+3, O_of_te             ,'-'  , ...
	     te*1E+3, gammaHz*Ge*ze_of_te ,'r--', ...
	     te*1E+3, DNu0_of_ze_of_te    ,'m-.');
	xlabel('t_e [ms]'); ylabel('O(t_e) [Hz]'); grid; set_gca;
    legend({'O(t_e)','\gamma*G_e*z_e(t_e)','\Delta\nu(z_e(t_e))'},'Location','Best');

    figure; hold; title('RF phase term due to \omega  vis-a-vis  the basic RF freqency');
	plot(te*1E+3, cumsum(2*pi*gammaHz*Ge*ze_of_te*dte) ,'-'  , ...
	     te*1E+3, cumsum(2*pi*DNu0_of_ze_of_te*dte)    ,'r--');
	xlabel('t_e [ms]'); ylabel('[rad]'); grid; set_gca;
    legend({'2*\pi*\int(gammaHz*Ge*ze_of_te)dt_e','2*\pi\int(DNu0_of_ze_of_te)dt_e'},'Location','Best');
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
ta  = Ta * (hi_res_za - Zia) / (Zfa - Zia);
dta = ta(2) - ta(1);
disp(sprintf('\n dte=%d\n dta=%d',dte,dta));
if (DEBUG_FLAG >= 2)
 	figure; hold; plot(hi_res_za,ta*1E+3,'-');
 	title('Acquisition temporal axis t_a(z)'); ylabel('t_a [ms]'); xlabel('z_a-axis [cm]'); grid; set_gca;
end;

% --------------------
% Acquisition gradient
% --------------------
% First calculate dDNu0 as a function of za.
% Since there's a linear (one-to-one) relation between za & ta it is sufficient to calculate d_DNu0/d_za
dDNu0_dza = polyval(P_dDNu0_dz,hi_res_za);
Ga_of_ta = -(1/gammaHz) * (2*alpha2*((Zfa-Zia)/Ta) + dDNu0_dza);
if (DEBUG_FLAG)
	figure; hold; plot(ta*1E+3,Ga_of_ta,'.-');
	title('Acquisition gradient G_a(t_a)'); xlabel('t_a-axis [ms]'); ylabel('G_a [G/cm]'); grid; set_gca;
end;

% ---------------------------------------------
% DEBUG: Total acquisition phase @ ta=0 & ta=Ta
% ---------------------------------------------
beta0  = ((Zfa-Zia)/Ta);
if (DEBUG_FLAG >= 3)
	plot_pre_acq_phase_1DUF180_RHR(z_axis_,Zia,alpha0,alpha1,alpha2,beta0,DNu0,INT_DNu0,OmegaE,Gpr,Ge,dDNu0_dz,Ga_of_ta,Phi_RF_of_te,context);
	plot_post_acq_phase_1DUF180_RHR(z_axis_,DNu0,OmegaE,Gpr,Ga_of_ta,context);
	plot_acq_phase_vs_t_1DUF180_RHR(z_axis_,DNu0,dDNu0_dz,OmegaE,Tp,te_of_z,Gpr,Tpr,ta,dta,Ga_of_ta);
end;

theoretical_phi_pre_acq1 = 2*pi*(alpha2*(z_axis.^2) + alpha1*z_axis + alpha0 - (1/beta0)*(DNu0.*(z_axis - Zia) - 2*INT_DNu0));
theoretical_phi_pre_acq2 = 2*pi*(alpha2*(z_axis.^2) + alpha1*z_axis + alpha0);
theoretical_phi_pre_acq3 = 2*pi*(- (1/beta0)*(DNu0.*(z_axis - Zia) - 2*INT_DNu0));

return;

