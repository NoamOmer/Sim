% Design the excitation pulse and acquisition gradient for fixing B0 inhomogeneity in 2D Hybrid UF/EPI
function [te,dte,O_of_te,Ge,Gpr,Ga_of_ta,R_of_te,chirp_rot,phi_chirp_rot] = design_inhom_fix_UF_1D180(context);
set_context;

z_axis_ = z_axis;   % Set the external 'z_axis' into the loacl 'z_axis_'

% ------------------
% Acquisition z-axis
% ------------------
% 'Z' axis for acquisition - These parems has to be defined here since we would like to be able to 
% use the current function without an external context
dz  = z_axis_(2)   - z_axis_(1);
Lz  = z_axis_(end) - z_axis_(1);
disp(sprintf('\n dz=%d;\n Lz=%d',dz,Lz));

Zia       = z_axis_(end);
Zfa       = z_axis_(1);
dZa       = dZa * (-sign(z_axis_(2) - z_axis_(1)));  % 'sign' is in case of reversed scan
hi_res_za = [Zia:(dZa/Za_AXIS_INTERP_FACTOR):Zfa];

if (DEBUG_FLAG)
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
if (DEBUG_FLAG)
	figure; hold; plot(z_axis_,DNu0,'b.-');
	title('Polynomial Approximation of Field Inhomogeneity \Delta\omega(z)');
	xlabel('z_e axis [cm]'); ylabel('\Delta\Omega [Hz]'); set_gca; grid;
end;

% ----------------
% Purge parameters
% ----------------
% Calculate alpha1 & 2 (alpha0 represents the RF initial phase and can be arbitrarily set? Maybe it relates to tof)
% Removed 2*pi from the SPA formula for the acquisition pixel size:
alpha2 = 1/(2*(dZa^2))                                 % [1/cm^2]
alpha1 = -(Ta/Lz)*(DNu0(end)) - 2*alpha2*z_axis_(end)  % [1/cm]   We assume that Zfe=Zia
alpha0 = 0;                                             % [none]

% ------------------------------
% Excitation duration & Gradient
% ------------------------------
mode = 3;
switch mode
    case 1  % Tp --> Ge
        Gpr = 5;
        Tp = (1/Ex_Fac)*( (Ta/Lz)*(DNu0(1)) + Ta*dDNu0_dz(1) + 2*alpha2*z_axis_(1) + alpha1 + 2*gammaHz*Gpr*Tpr ) /...
            (dDNu0_dz(end) - dDNu0_dz(1));
        Ge = -(1/(gammaHz*Tp)) * ...
            ((Ta/Lz)*DNu0(1) + (Tp+Ta)*dDNu0_dz(1) + 2*alpha2*z_axis_(1) + alpha1 + gammaHz*Gpr*Tpr);
    case 2  % Ge --> Tp
        % Exactly the same solution of 2 equations with 2 unknowns - except - solve Ge 1st and Tp 2nd
        % Gives the same results
        Gpr = 5;
        LHS = 2*alpha2*z_axis_(1) + alpha1 + (Ta + Tpr)*dDNu0_dz(1) + (Ta/Lz)*DNu0(1);
        Ge = (-LHS*dDNu0_dz(end) - (gammaHz*Gpr + dDNu0_dz(end))*Tpr*dDNu0_dz(1)) / ...
             ((gammaHz)*(LHS + (gammaHz*Gpr + dDNu0_dz(end))*Tpr));
        Tp = (gammaHz*Gpr + dDNu0_dz(end))*Tpr / (gammaHz*Ge + dDNu0_dz(end));
    case 3  % Ge --> Gpr
        Gpr = (1/(2*gammaHz*Tpr)) * ...
              (Tp*(dDNu0_dz(end)-dDNu0_dz(1)) - (Ta/Lz)*DNu0(1) - Ta*dDNu0_dz(1) - 2*alpha2*z_axis_(1) - alpha1);
        Ge  = -(1/(gammaHz*Tp))   * ...
              ((Ta/Lz)*DNu0(1) + (Tp+Ta)*dDNu0_dz(1) + 2*alpha2*z_axis_(1) + alpha1 + gammaHz*Gpr*Tpr);
end;
disp(sprintf('\n Tpr=%d\n Gpr=%d\n Tp=%d\n Ge=%d\n Ta=%d\n',Tpr,Gpr,Tp,Ge,Ta));

if (Inhomo_Mult_Factor == 0)
    if (abs(Ge - alpha2*Lz/(gammaHz*Tp)) > abs(Ge/10))
        error('Ge is not equal to alpha2*Lz/(gammaHz*Tp) for zero inhomogeneity');
    end;
end;

% --------------------
% Excitation time axis
% --------------------
% Calculate the effective frequency of spins as a fuction of 'z' during excitation: Omega(z)
OmegaE = gammaHz*Ge*z_axis_ + DNu0 + omega_cs;
% >>> VERIFICATION ------------------ %
OmegaE_min = min(OmegaE);
OmegaE_max = max(OmegaE);
OmegaE_axis  = linspace(OmegaE_min,OmegaE_max,length(z_axis_));
ze_of_OmegaE = interp1(OmegaE,z_axis_,OmegaE_axis);
dze_dOmegaE_app = diff(ze_of_OmegaE)/(OmegaE_axis(2)-OmegaE_axis(1));  % Approximated value
% ------------------------------- <<< %
dOmegaE_dz = gammaHz*Ge + dDNu0_dz;
dz_dOmegaE = 1./dOmegaE_dz;
if (DEBUG_FLAG >= 2)
	figure; hold on;
	subplot(2,3,1); plot(z_axis_,OmegaE,'.-',z_axis_,gammaHz*Ge*z_axis_,'r.-');
    xlabel('z-axis [cm]');   ylabel('[Hz]');      title('\omega_E(z)');    grid;   set_gca;
	legend({'\omega_0+\Delta\omega_0','\omega_0'},'Location','Best');
    
    subplot(2,3,2); plot(z_axis_,dOmegaE_dz,'.-');
    xlabel('z-axis [cm]');   ylabel('[Hz/cm]');   title('d\omega_E/dz');   grid;   set_gca;

    subplot(2,3,3); plot(z_axis_,dz_dOmegaE,'.-');
    xlabel('z-axis [cm]');   ylabel('[cm/Hz]');   title('dz/d\omega_E');   grid;   set_gca;

    % >>> VERIFICATION ------------------ %
	subplot(2,3,4); plot(OmegaE_axis,ze_of_OmegaE,'.-',gammaHz*Ge*z_axis_,z_axis_,'r.-'); xlabel('[Hz]'); ylabel('[cm]');
	legend({'z(\omega_0+\Delta\omega_0)','z(\omega_0)'},'Location','Best'); title('z(\omega_E) (app)');   grid;    set_gca;
	subplot(2,3,5); plot(ze_of_OmegaE(2:end),1./dze_dOmegaE_app,'.-'); xlabel('z-axis [cm]');  ylabel('[Hz/cm]');
	title('d\omega_E/dz (app)');      grid;    set_gca;
	subplot(2,3,6); plot(ze_of_OmegaE(2:end),dze_dOmegaE_app,'.-');    xlabel('z-axis [cm]');  ylabel('[cm/Hz]');
	title('1/(d\omega_E/dz) (app)');  grid;    set_gca;
	% --------------------------------<<< %
end;
% Calculate the relation between time and location during excitation: te(z)
te_of_z = 0.5 * Tp + ...
          0.5 * dz_dOmegaE .* ( 2*alpha2*z_axis_ + alpha1 - ...
                                (Ta/Lz) * (dDNu0_dz.*(z_axis_ - z_axis_(end)) - DNu0) + ...
                                gammaHz*Gpr*Tpr );
% Calculate a linear te axis
te  = linspace(min(te_of_z),max(te_of_z),length(z_axis_)); % Verify: should be 0 --> Tp
dte = te(2) - te(1);

if (DEBUG_FLAG)
	figure; hold; title('Excitation time of z-slices (t_e(z))');
	plot(z_axis_,te_of_z*1E+3,'.-',z_axis_,te*1E+3,'g+');
	xlabel('z-axis [cm]'); ylabel('t_e [ms]'); grid; set_gca;
	legend({'t_e(z)','t_e axis'},'Location','Best');
end;

% ----------------
% Excitation phase
% ----------------
% Calculate ze(te) during excitation by reversing te(z)
ze_of_te = interp1(te_of_z,z_axis_,te);  % 'ze' at equaly spaces points along time axis 'te'.
if (DEBUG_FLAG)
	figure; hold; title('z-Slice location along excitation time axis (z_e(t_e))');
	plot(te*1E+3,ze_of_te,'.-',te*1E+3,z_axis_,'g+');
	xlabel('t_e [ms]'); ylabel('z_e-axis [cm]'); grid; set_gca;
	legend({'z_e(t_e)','z_e axis'},'Location','Best');
end;
% sz_factor = round(0.075*length(ze_of_te));
% dze_of_te_start = ze_of_te(2) - ze_of_te(1);
% ze_of_te_start  = ze_of_te(1) - sz_factor*dze_of_te_start;
% dze_of_te_end   = ze_of_te(end) - ze_of_te(end-1);
% ze_of_te_end    = ze_of_te(end) + sz_factor*dze_of_te_end;
% ze_of_te = [ze_of_te_start:dze_of_te_start:ze_of_te(1),ze_of_te(2:end-1),ze_of_te(end):dze_of_te_end:ze_of_te_end];
% if (DEBUG_FLAG)
% 	figure; hold; title('z-Slice location along excitation time axis (z_e(t_e)) - Extrapolated');
% 	plot(te*1E+3,ze_of_te,'.-',te*1E+3,z_axis_,'g+');
% 	xlabel('t_e [ms]'); ylabel('z_e-axis [cm]'); grid; set_gca;
% 	legend({'z_e(t_e)','z_e axis'},'Location','Best');
% end;

% Calculate O(te) - the temporal derivative of the RF phase
% Note that we first need to calculate DNu0 as a function of ze (as opposed to DNu0(z))
DNu0_of_ze_of_te = polyval(P_DNu0,ze_of_te);
O_of_te = -(gammaHz*Ge*ze_of_te + DNu0_of_ze_of_te + omega_cs);        % [Hz]
Phi_RF_of_te = 2*pi*dte*cumsum(O_of_te);                               % [rad]
% Suggestion: maybe interpolate O(te) in order to get a more accurate result of the integral
if (DEBUG_FLAG)
	figure; hold; title('d\Phi_R_F(t_e)/dt = O(t_e)');
	plot(te*1E+3,O_of_te,'-'); xlabel('t_e [ms]'); ylabel('O(t_e) [Hz]'); grid; set_gca;
	plot(te*1E+3,-gammaHz*Ge*ze_of_te,'r--'); xlabel('t_e [ms]'); ylabel('O(t_e) [Hz]');
	plot(te*1E+3,-DNu0_of_ze_of_te   ,'m-.'); xlabel('t_e [ms]'); ylabel('O(t_e) [Hz]');
    legend({'O(t_e)','-\gamma*G_e*z_e(t_e)','-\Delta\nu(z_e(t_e))'},'Location','Best');
end;

% ----------------
% Excitation power
% ----------------
R_of_te = numeric_derivation(O_of_te,te);    % [1/sec^2]    Chirp rate (dO/dte)
% Note that the following function takes omega_cs into consideration. This means that it is accounted
% for twice instead of once, since it is also explicitly added to O_of_te.
[chirp_rot, phi_chirp_rot] = gen_amp_n_phase_modulated_Chirp_pulse(te,R_of_te,Phi_RF_of_te,context); % [[T],[rad]]
if (DEBUG_FLAG)
	figure; hold;
	plot(te*1E+3,phi_chirp_rot/2*pi,'.-'); title('\Phi_R_F (t_e) [unitless]'); xlabel('Time [ms]'); ylabel('[none]');
	plot(te*1E+3,unwrap(angle(chirp_rot))/2*pi,'g.-'); legend({'\phi_R_F','unwrap(angle(RF))'},'Location','Best');
end;

% ---------------------
% Acquisition time axis
% ---------------------
% Create a time axis for acquisition
ta  = -(hi_res_za - hi_res_za(1)) * (Ta/Lz);
dta = ta(2) - ta(1);
disp(sprintf('\n dte=%d\n dta=%d',dte,dta));
if (DEBUG_FLAG)
 	figure; hold; plot(hi_res_za,ta*1E+3,'.-');
 	title('Acquisition temporal axis t_a(z)'); ylabel('t_a [ms]'); xlabel('z_a-axis [cm]'); grid; set_gca;
end;

% --------------------
% Acquisition gradient
% --------------------
% First calculate dDNu0 as a function of za.
% Since there's a linear (one-to-one) relation between za & ta it is sufficient to calculate d_DNu0/d_za
dDNu0_dza = polyval(P_dDNu0_dz,hi_res_za);
Ga_of_ta = -(1/gammaHz) * (2*alpha2*(Lz/Ta) + dDNu0_dza);
if (DEBUG_FLAG)
	figure; hold; plot(ta*1E+3,Ga_of_ta,'.-');
	title('Acquisition gradient G_a(t_a)'); xlabel('t_a-axis [ms]'); ylabel('G_a [G/cm]'); grid; set_gca;
end;

% ---------------------------------------------
% DEBUG: Total acquisition phase @ ta=0 & ta=Ta
% ---------------------------------------------
if (DEBUG_FLAG >= 2)
	plot_pre_acq_phase_1DUF180(z_axis_,alpha0,alpha1,alpha2,DNu0,INT_DNu0,OmegaE,Gpr,Ge,dDNu0_dz,Ga_of_ta,Phi_RF_of_te,context);
	plot_post_acq_phase_1DUF180(z_axis_,DNu0,OmegaE,Gpr,Ga_of_ta,context);
	plot_acq_phase_vs_t_1DUF180(z_axis_,DNu0,OmegaE,Tp,te_of_z,Gpr,Tpr,ta,dta,Ga_of_ta);
end;

return;

