% Design the excitation pulse and acquisition gradient for fixing B0 inhomogeneity in 1DUF
% Uses Right-Hand Spin Rotation
% Ver_7  : LHR --> RHR
% Ver_8  : Discarded
% Ver_9  : Instead of determining dZa, determine Ge
% Ver_9.1: Adjust the sequence for a 2D Hybrid imaging:
%           - We still assume that the acquisition is continuous
%           - We add into consideration the 180 pulse w/o a balancing delay
%           - We divide the purge stage into two stages pre- & post-180 and thus
%             do not take it into consideration in the formalism
%           - We still assume (in the last stages of the algorithm) that Zia=Zfe & Zfa=Zie

%  Output parameters:
%  ------------------
%  te            :  1D vector            :  [sec]    :   Time axis
%  O_of_te       :  1D vector (te)       :  [Hz]     :   RF pulse time-dependant phase (dPhiRF_dte)
%  Ge            :     scalar            :  [G/cm]   :   Excitation gradient
%  ta            :  1D vector (za_axis)  :  [sec]    :   Acquisition temporal axis
%  Ga_of_ta      :  1D vector (za_axis)  :  [G/cm]   :   Acquisition gradient
%  R_of_te       :  1D vector            :  [Hz/sec] :   Time dependent chirp rate
%  chirp_rot     :                       :  []       :   
%  phi_chirp_rot :                       :  []       :   
%
%  Internal parameters:
%  --------------------
%  z_axis_       :  1D vector            :  [cm]     :   Sample spatial axis

function [te_,dte_,O_of_te,Ge,ta,Ga_of_ta,R_of_te,chirp_rot,phi_chirp_rot] = design_inhom_fix_1DUF_RHR(context);
set_context;
% Set the external 'z_axis' into the loacl 'z_axis_'
z_axis_ = z_axis;

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
hi_res_za = linspace(Zia,Zfa,n_Ga_pts);

if (DEBUG_FLAG >= 2)
	figure; hold; plot(1:length(z_axis_),z_axis_,'b.-');
	title('Excitation spatial axes');  xlabel('[none]');  ylabel('[cm]');
	legend({'z_e axis'},'Location','Best'); set_gca; grid;
end;

% ----------------------------------------------------------------
% Calculate the needed derivations (and integrals) of Delta_Nu0(z)
% ----------------------------------------------------------------
P_INT_DNu0 = polyint(P_DNu0);
P_dDNu0_dz = polyder(P_DNu0);

DNu0     = polyval(P_DNu0     ,z_axis_);     % [Hz]     Polynomial approximation for B0 inhomogeneity
INT_DNu0 = polyval(P_INT_DNu0 ,z_axis_);     % [Hz*cm]  spatial integral
dDNu0_dz = polyval(P_dDNu0_dz ,z_axis_);     % [Hz/cm]  1st spatial derivative
if (DEBUG_FLAG >= 1)
	fh1 = figure; hold on; plot(z_axis_,DNu0,'b-');
	title('Inhomogeneity Poly-fit \Delta\omega (z)');
	xlabel('z_e axis [cm]'); ylabel('\Delta\Omega [Hz]'); set_gca; grid;
% 	fh2 = figure; hold on; plot(z_axis_,INT_DNu0,'b-');
% 	title('Integral of Inhomogeneity Poly-fit \int\Delta\omega (z)dz');
% 	xlabel('z_e axis [cm]'); ylabel('\int\Delta\Omega(z)dz [Hz*cm]'); set_gca; grid;
% 	fh3 = figure; hold on; plot(z_axis_,dDNu0_dz,'b-');
% 	title('Derivative of Inhomogeneity Poly-fit d\Delta\omega(z)/dz');
% 	xlabel('z_e axis [cm]'); ylabel('d\Delta\Omega/dz [Hz/cm]'); set_gca; grid;
end;
DNu0     = smooth_edges(DNu0    ,context);
INT_DNu0 = smooth_edges(INT_DNu0,context);
dDNu0_dz = smooth_edges(dDNu0_dz,context);
if (DEBUG_FLAG >= 2)
	figure(fh1); hold on; plot(z_axis_,DNu0,'r-');      legend({'pre-smooth','post-smooth'},'Location','Best');
% 	figure(fh2); hold on; plot(z_axis_,INT_DNu0,'r-');  legend({'pre-smooth','post-smooth'},'Location','Best');
% 	figure(fh3); hold on; plot(z_axis_,dDNu0_dz,'r-');  legend({'pre-smooth','post-smooth'},'Location','Best');
end;

% ---------------------------
% Quadratic phase coefficient
% ---------------------------
alpha2 = (1/(2*(Zie_-Zia))) * ( (Ta/(Zfa-Zia)) * (DNu0(end) - DNu0(1)) + ...
                                (Ta-Tp)*dDNu0_dz(1)                    - ...
                                gammaHz*Ge*Tp );                             % [1/cm^2]
alpha1 = -(Ta/(Zfa-Zia))*(DNu0(end)) - 2*alpha2*Zia;                         % [1/cm]
alpha0 = 0;                                                                  % [none]
dZa = sign(alpha2) * 1/sqrt(2*abs(alpha2));            % [cm] w/o inhomogeneity = sqrt(L/(gammaHz*Ge*Tp))

disp(sprintf('dZa = %5.3f [mm]'  ,dZa*10));
disp(sprintf('Ge  = %5.3f [G/cm]',Ge    ));

% --------------------
% Excitation time axis
% --------------------
% Calculate the effective frequency of spins as a fuction of 'z' during excitation: Omega(z)
OmegaE     = gammaHz*Ge*z_axis_ + DNu0 + omega_cs;
dOmegaE_dz = gammaHz*Ge         + dDNu0_dz;
dz_dOmegaE = 1./dOmegaE_dz;

OmegaE_axis  = linspace(min(OmegaE),max(OmegaE),length(z_axis_));
ze_of_OmegaE = interp1(OmegaE,z_axis_,OmegaE_axis);                 % Required for debugging and for SLR
if (DEBUG_FLAG >= 2)
	dze_dOmegaE_app = numeric_derivation(ze_of_OmegaE,OmegaE_axis); % Approximated value
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
te_of_z = ((2*alpha2*z_axis_ + alpha1 - (Ta/(Zfa-Zia))*(dDNu0_dz.*(z_axis_ - Zia) - DNu0)) ./ dOmegaE_dz) + Tp;
if (DEBUG_FLAG >= 2)
	figure; hold; title('Excitation time of z-slices (t_e(z))');
	plot(z_axis_,te_of_z*1E+3,'-'); plot([z_axis_(1) z_axis_(end)], [te_of_z(1) te_of_z(end)]*1E+3, 'k--');
	xlabel('z-axis [cm]'); ylabel('t_e [ms]'); grid; set_gca;
end;
if ((sum(te_of_z < 0) > 1)         || (min(te_of_z) ~= te_of_z(1)) || ...
    (max(te_of_z) ~= te_of_z(end)) || (~isempty(find(diff(te_of_z) <= 0))))
	figure; hold; title('Excitation time of z-slices (t_e(z))');
	plot(z_axis_,te_of_z*1E+3,'-'); plot([z_axis_(1) z_axis_(end)], [te_of_z(1) te_of_z(end)]*1E+3, 'k--');
	xlabel('z-axis [cm]'); ylabel('t_e [ms]'); grid; set_gca;
	h = msgbox('te(z) has illegal values. Would you like to continue?');
	uiwait(h);
end;

% Calculate a linear te axis
te_  = linspace(te_of_z(1),te_of_z(end),length(z_axis_));
dte_ = te_(2) - te_(1);

% ----------------
% Excitation phase
% ----------------
% Calculate ze(te) during excitation by reversing te(z)
ze_of_te = interp1(te_of_z,z_axis_,te_);  % 'ze' at equaly spaces points along time axis 'te'.
if (DEBUG_FLAG >= 2)
	figure; hold; title('z-Slice location along excitation time axis (z_e(t_e))');
	plot(te_*1E+3,ze_of_te,'-'); plot([te_(1) te_(end)]*1E+3, [ze_of_te(1) ze_of_te(end)], 'k--');
	xlabel('t_e [ms]'); ylabel('z_e-axis [cm]'); grid; set_gca;
end;

% Calculate the RF phase temporal derivation O(te).
% Note that we first need to calculate DNu0 as a function of ze (as opposed to DNu0(z))
DNu0_of_ze_of_te = polyval(P_DNu0,ze_of_te);                     % [Hz]
% figure; hold on; plot(ze_of_te,DNu0_of_ze_of_te,'b-');
% DNu0_of_ze_of_te = smooth_edges(DNu0_of_ze_of_te,context);
% hold on; plot(ze_of_te,DNu0_of_ze_of_te,'r-'); title('DNu0_of_ze_of_te');

O_of_te = (gammaHz*Ge*ze_of_te + DNu0_of_ze_of_te + omega_cs);   % [Hz]  Chirp frequency (te)
Phi_RF_of_te = cumsum(2*pi*O_of_te*dte_);                         % [rad]
% Suggestion: maybe interpolate O(te) in order to get a more accurate result of the integral
if (DEBUG_FLAG >= 1)
    f1 = 1e-3*gammaHz*Ge*ze_of_te;
    f2 = 1e-3*DNu0_of_ze_of_te;
    figure;
    subplot(2,2,1); plot(te_*1E+3,f1   ,'.-'); title('gammaHz*Ge*ze(te)'); xlabel('t_e [ms]'); ylabel('[kHz]'); grid; set_gca;
    subplot(2,2,2); plot(te_*1E+3,f2   ,'.-'); title('DNu0(ze(te))'     ); xlabel('t_e [ms]'); ylabel('[kHz]'); grid; set_gca;

% 	figure; hold; title('O(t_e) = d\Phi_R_F/dt_e: Contribution of G_e vis-a-vis \Delta\nu');
% 	plot(te_*1E+3, 1E-3*O_of_te             ,'-'  , ...
% 	     te_*1E+3, 1E-3*gammaHz*Ge*ze_of_te ,'r--', ...
% 	     te_*1E+3, 1E-3*DNu0_of_ze_of_te    ,'m-.');
% 	xlabel('t_e [ms]'); ylabel('O(t_e) [kHz]'); grid; set_gca;
%     legend({'O(t_e)','\gamma*G_e*z_e(t_e)','\Delta\nu(z_e(t_e))'},'Location','Best');
% 
%     figure; hold; title('RF phase term due to \omega  vis-a-vis  the basic RF freqency');
% 	plot(te_*1E+3, cumsum(2*pi*gammaHz*Ge*ze_of_te*dte_) ,'-'  , ...
% 	     te_*1E+3, cumsum(2*pi*DNu0_of_ze_of_te*dte_)    ,'r--');
% 	xlabel('t_e [ms]'); ylabel('[rad]'); grid; set_gca;
%     legend({'2*\pi*\int(\gamma*Ge*z_e(t_e))dt_e','2*\pi\int(DNu0(z_e(t_e)))dt_e'},'Location','Best');
end;

% ----------------
% Excitation power
% ----------------
R_of_te = numeric_derivation(O_of_te,te_);    % [1/sec^2]    Chirp rate (dO/dte)
% !!!  Note  !!!
% The following API takes omega_cs into consideration. This means that it is accounted
% for twice instead of once, since it is also explicitly added to O_of_te.
[chirp_rot, phi_chirp_rot] = gen_amp_n_phase_modulated_Chirp_pulse(te_,R_of_te,Phi_RF_of_te,context); % [[T],[rad]]

% ---------------------
% Acquisition time axis
% ---------------------
% Create a time axis for acquisition
ta  = Ta * (hi_res_za - Zia) / (Zfa - Zia);
dta = ta(2) - ta(1);
if (DEBUG_FLAG >= 2)
	figure; hold; plot(ta*1E+3,hi_res_za,'-');
	title('Acquisition spatial axis z_a(t_a)'); xlabel('t_a [ms]'); ylabel('z_a-axis [cm]'); grid; set_gca;
end;

% -------------------- %
% Acquisition gradient %
% -------------------- %
% First calculate dDNu0 as a function of za.
% Since there's a linear (one-to-one) relation between za & ta it is sufficient to calculate d_DNu0/d_za
dDNu0_dza = polyval(P_dDNu0_dz,hi_res_za);
% figure; hold on; plot(hi_res_za,dDNu0_dza,'b-');
% dDNu0_dza = smooth_edges(dDNu0_dza,context);
% hold on; plot(hi_res_za,dDNu0_dza,'r-'); title('dDNu0_dza');

Ga_of_ta = -(1/gammaHz) * (2*alpha2*((Zfa-Zia)/Ta) + dDNu0_dza);
if (DEBUG_FLAG >= 2)
	figure; hold; plot(ta*1E+3,Ga_of_ta,'-');
	title('Acquisition gradient G_a(t_a)'); xlabel('t_a-axis [ms]'); ylabel('G_a [G/cm]'); grid; set_gca;
    
    g1 = -(1/gammaHz) * (2*alpha2*((Zfa-Zia)/Ta));
    g2 = -(1/gammaHz) * (dDNu0_dza);
    figure;
    subplot(2,1,1); plot(        g1,'.-'); title('g1'); xlabel('t_a-axis [ms]'); ylabel('[G/cm]'); grid; set_gca;
    subplot(2,1,2); plot(ta*1E+3,g2,'.-'); title('g2'); xlabel('t_a-axis [ms]'); ylabel('[G/cm]'); grid; set_gca;
%    figure; plot(hi_res_za,dDNu0_dza,'.-',z_axis_,dDNu0_dz,'.-'); legend({'dDNu0_dza','dDNu0_dz'});
end;

% ---------------------------------------------
% DEBUG: Total acquisition phase @ ta=0 & ta=Ta
% ---------------------------------------------
if (DEBUG_FLAG >= 5)
	pack; plot_pre_acq_phase_1DUF_RHR(z_axis_,Zia,Zfa,alpha0,alpha1,alpha2,DNu0,dDNu0_dz,INT_DNu0,OmegaE,te_of_z,Ge,O_of_te,Phi_RF_of_te,Tp,Ta);
	pack; plot_post_acq_phase_1DUF_RHR(z_axis_,DNu0,OmegaE,te_of_z,Ga_of_ta,context);
	pack; plot_acq_phase_vs_t_1DUF_RHR(z_axis_,DNu0,dDNu0_dz,OmegaE,Tp,te_of_z,ta,dta,Ga_of_ta,Ge);
end;

% -----------------------------------------------------------
% Ugly: necessary for external (mainly verification) purpuses
% -----------------------------------------------------------
% 'te' as a function of EXTERNAL 'z' axis -
% te_of_z_axis = interp1(z_axis_  ,te_of_z,z_axis);    % <-- uses EXTERNAL z_axis

% 'ta' as a function of EXTERNAL 'z' axis -
ta_of_z = interp1(hi_res_za, ta, z_axis);              % <-- uses EXTERNAL z_axis

% Required post-excitation phase profile
Phi_e   = alpha2*(z_axis_.^2) + alpha1*z_axis_ + alpha0 - (Ta/(Zfa-Zia)) * (DNu0 .*(z_axis_ - Zia) - 2*INT_DNu0);

% Quadratic part of post-excitation phase profile - only in perfectly homogeneous field
Phi_e_sq = -alpha2*(z_axis_.^2);
% dbstop in design_inhom_fix_1DUF_RHR at 252

return;

function [smoothed_vec] = smooth_edges(vec,context)
set_context;

initial_vec = vec;
dvec = numeric_derivation(vec,1:length(vec));
tmp1 = length(dvec);
tmp2 = [(-round(tmp1/2)+1):1:(+round(tmp1/2))];
tmp2 = tmp2(1:tmp1);
dvec_win = exp(-(2.25*abs(tmp2./tmp1)).^10);
% figure;plot(dvec_win);
dvec = dvec .* dvec_win;
int_dvec = cumsum(dvec);
smoothed_vec = int_dvec;
smoothed_vec = smoothed_vec + vec(round(length(vec)/2)) - smoothed_vec(round(length(smoothed_vec)/2));

% smoothed_vec = transpose(smooth(smoothed_vec,20,'loess'));
smoothed_vec = initial_vec; % do nothing

% ---- %
% JUNK %
% ---- %
% initial_vec = vec;
% % dbstop in design_inhom_fix_1DUF_RHR.m at 265
% p = 7.5; % actually, double that precentage
% N = length(vec);
% n = round(length(vec)/p);
% win  = round(n/2);
% idx1 = 1; % n - round(1.7*win);
% idx2 = n + round(1.0*win);
% idx3 = N - n + 1 - round(1.0*win);
% idx4 = N; % N - n + 1 + round(1.7*win);
% vec(1:n)          = vec(n);        vec1 = smooth(vec,win);  vec(idx1:idx2) = vec1(idx1:idx2);
% vec(N:-1:(N-n+1)) = vec((N-n+1));  vec1 = smooth(vec,win);  vec(idx3:idx4) = vec1(idx3:idx4);

% smoothed_vec = initial_vec; % do nothing
% smoothed_vec = vec;
return;


