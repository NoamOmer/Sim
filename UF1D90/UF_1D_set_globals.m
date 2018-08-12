
set_globals;
global sim_name;         sim_name         = 'UF 1D Simulation';
global Seq1D;            Seq1D            = 1; 
global Seq2D;            Seq2D            = 0;

global Inhomo_Flag;      Inhomo_Flag      =	1;                                % A flag: enable/disable B0 inhomogeneity
global DesignPulse_Flag; DesignPulse_Flag = 0;                                % A flag: enable/disabeB0 inhomogeneity fix using designed pulses
global Inhomo_Fix_Flag;  Inhomo_Fix_Flag  = 0;                                % A flag: enable/disabeB0 inhomogeneity fix
global RFnGRAD_Win_Flag; RFnGRAD_Win_Flag = 1;                                % A flag: enable/disabe designed Chirp windowing
global SLR_flag;         SLR_flag         = 0;                                % A flag: enable/disabe SLR pulse
global SE_flag;          SE_flag          = 1;                                % A flag: enable/disabe SE PI pulse

global M0;               M0               = [0,0,1];                          % [amper/m] Magnetization
global M_init;
global B0_lab;           B0_lab           = [0,0,3];                          % [Tesla]   Lab magnetic field

global RF_ESP;           RF_ESP           = 500;
global RF_ESP_SLR;       RF_ESP_SLR       = 500;

%  -----------------------------
%   Number of stationary points
%  -----------------------------
global Nsp;              Nsp              = 1;
global FIDnum;
if (Nsp > 1) && (isempty(FIDnum))
	FIDnum = 600;
	fn = sprintf('D:\\PhD\\Matlab\\Simulations\\Ver_4\\UF1D90\\MultSP_sim\\%3.0f.fid',FIDnum);
	while (exist(fn,'file') == 2)
		op.Resize      = 'on';
		op.WindowStyle = 'normal';
		FIDnumCell = inputdlg('Please enter a new FID number',sprintf('File %3.0f exist',FIDnum),1,{'FID #'},op);
		if (isempty(FIDnumCell))
			FIDnum = FIDnum + 1;
			uiwait(msgbox(sprintf('Incrementing FIDnum by 1 to %3.0f',FIDnum)));
		else
			FIDnum = str2double(FIDnumCell{1});
		end;
		fn = sprintf('D:\\PhD\\Matlab\\Simulations\\Ver_4\\UF1D90\\MultSP_sim\\%3.0f.fid',FIDnum);
	end;
end;

%  ===================
%   Manual parameters
%  ===================
global sample_Nz;        sample_Nz        = 15e4;
global required_res;     required_res     = 0.20;                               % [cm]
global BWppacq;          BWppacq          = 333;                                % [Hz/Px] Acquisition duration
global Ta;               Ta               = (1/BWppacq)/Nsp;                    % [sec]   Acquisition duration
global ROOS;             ROOS             = 1.0;
global Lhalf;            Lhalf            = 2.56*6;                             % [cm]
global Tcrush;           Tcrush           = 2000e-6;                            % [sec]   Crusher gradient duration
global Gcrush;           Gcrush           = 0.00;                               % [G/cm]  Crusher gradient power

global exc_phase;        exc_phase        = pi/2;
% dte_interp_F = 2;                                                             %         Just to be on the safe side, divide required dte by 2

global Trefoc;           Trefoc           = 3000e-6;                            % [sec]
global N_refoc;          N_refoc          = 600;
global n_lobes_refoc;    n_lobes_refoc    = 2;
global Grefoc;           Grefoc           = 0;                                  % [G/cm]  Refocusing gradient
global refoc_angle;      refoc_angle      = 180;                                % [deg]   Refocusing flip angle
global refoc_phase;      refoc_phase      = 0;                                  % [rad]   Refocusing phase


global Tp;               Tp               = Ta;                                 % [sec]   Duration of Chirp pulse
global Dk;               Dk               = 1/required_res;                     % [cm]
global Ga;               Ga               = [0,0,Dk / (gammaHz*Ta)];            % [G/cm]  Acquisition gradient
global Ge;               Ge               = Ga;                                 % [G/cm]  Magnetic field gradient
global fb;               fb               = abs(gammaHz*Ga(3)*Lhalf*2);         % [Hz]    Filter bandwidth
global dta;              dta              = 1.0*1/fb;                           % [sec]   Acquisition dwell time
global N_acq_pts;        N_acq_pts        = ceil(Ta/dta);

%  ---------------------
%   'Z' EXCITATION axis
%  ---------------------
global Zie;              Zie              = -Lhalf;                             % [cm] Excitation z-axis beginning
global Zfe;              Zfe              = +Lhalf;                             % [cm] Excitation z-axis end
global Zia;              Zia              = Zfe;
global Zfa;              Zfa              = Zie;
global Lz_imag;          Lz_imag          = abs(Zfe-Zie);                       % [cm] Excitation length
global z_axis;                                                                  % [cm] EXCITATION (!) z-axis 
global alpha1;
global alpha2;


%  -------------
%   Sample axis
%  -------------
Shalf = Lhalf;
global sample_Zi;        sample_Zi        = +Shalf*sign(Zie);                  % [cm]      Sample z-axis
global sample_Zf;        sample_Zf        = +Shalf*sign(Zfe);
global sample_Lz;        sample_Lz        = abs(sample_Zi-sample_Zf);
global sample_dz;        sample_dz        = sample_Lz / sample_Nz;
global sample_z_axis;    sample_z_axis    = linspace(sample_Zi, sample_Zf, sample_Nz);


%  ---------------
%   Inhomogeneity
%  ---------------
global Map_Mult_Factor;  Map_Mult_Factor  = 0;                                % A flag: enable/disabe zero inhomogeneity test
global gPolyOrder;       gPolyOrder       = 3;
global P_DNu0;
global dOmega_0;                                                              % [Hz] dOmega_0 = Eta0+Eta1*z_axis+Eta2*(z_axis.^2)+Eta3*(z_axis.^3);
global dB0z;                                                                  % [T]

%  ------------
%   Excitation
%  ------------
global rfwdth;
if (SLR_flag)
	rfwdth = Tp + 2E-3;                                                       % [sec]
else
	rfwdth = Tp + 0E-3;                                                       % [sec]
end;
global dte;                                                                   % [sec]     Temporal resolution
global te;                                                                    % [sec]     Time axis
global te_of_z;                                                               % [sec]     Time axis vs. spatial location 'z'
global O_of_te;                                                               % [Hz]      O(kHz)
global Oi_Hz;                                                                 % [Hz]      O(kHz)
global rate;                                                                  % [1/sec^2] O(kHz/ms)

global OmegaE;
global ze_of_OmegaE;
global Phi_e;                                                                 % [none]    Phase at end of excitation

% B_RF amplitude in rot. frame: on the 300MHz the maximal value is 8kHz
% global B1MaxampLabHz;  B1MaxampLabHz = 16E+3;                               % [Hz]      RF field lab maximal amplitude
global B1MaxampLabHz;    B1MaxampLabHz    = 2*22e-6*gamma_Hz_T;               % [Hz]      Siemens
global B1ampLabHz;                                                            % [Hz]      Non time-modulated RF amp O(kHz)
global B1ampCoeff;                                                            % [none]    B1 amp coefficient

% 0.26 constant is given in JMR 172 (2005) 179–190.
% We use 0.52 for the lab frame and then divide it by two for the rotating frame.
B1ampCoeff = 0.265; % Determined manually (just run the simulation and check at what value the excitation is exactly pi/2)
if (DesignPulse_Flag)
	% Ge, te & O are set by inhomogeneity correction algorithm
else
	Oi_Hz  = gammaHz*Ge(3)*Zie;
	Of_Hz  = gammaHz*Ge(3)*Zfe;
	DO_Hz  = (Of_Hz - Oi_Hz)/Nsp;
	rate   = DO_Hz/Tp;
	DO_Hz_abs = abs(DO_Hz);
% 	DO_kHz = 1E-3*DO_Hz;
	dte = 1/(2*DO_Hz_abs);
% 	dte = dte/dte_interp_F;
	two_powers=[128,256,512,1024,2048,3082,4096,8192];
	n_RF_pts = Tp/dte;
	loc = find(two_powers > n_RF_pts);
	n_RF_pts = two_powers(loc(1));
	dte = Tp/n_RF_pts;
	
	te  = 0:dte:(rfwdth-dte);
	B1ampLabHz = B1ampCoeff*sqrt(abs(rate));
end;


%  -------------
%   Acquisition
%  -------------
% 'Z' acquisition axis
global dZa;                                                                   % [cm]   Effective acquisition resolution
global za_axis;                                                               % [cm]   Acquisition z-axis
global hi_res_za;                                                             % [cm]   High res. acq. z-axis (matching 'ta')
global ta;                                                                    % [sec]  Time axis vs. hi_res_za    (Zia-->Zfa)
global ta_of_z;                                                               % [sec]  Time axis vs. basic z_axis (Zie-->Zfe)
% global GRADdt;           GRADdt           = 4e-6;                             % [sec]  Acquisition temporal resolution
global Ga_of_ta;                                                              % [G/cm] Time-dependent acquisition gradient
if (DesignPulse_Flag)
	% ta, dta & Ga are set by inhomogeneity correction algorithm
else
	ta  = 0:dta:(Ta-dta);
% 	Ga  = (Ge*Tp/Ta) * sign(SE_flag-0.5);    % reverse Ga in case a SE was performed
	Ga  = Ga * sign(SE_flag-0.5);            % reverse Ga in case a SE was performed
	dZa = sqrt((Lz_imag/Nsp)/(gammaHz*abs(Ge(3))*Tp));
	dZa = sqrt(rate)/(gammaHz*abs(Ge(3)));
end;
global F_SR; F_SR = N_acq_pts / (Lz_imag/dZa);

% Verify simulation spatial axis
if (2*pi*sample_dz*gammaHz*Ge(3) > pi/1)  % pi/3 instead of pi or 2*pi just to be on the safe side.
	error(sprintf('sample spatial resolution is too low (%3.1f > %3.1f) [Ge=%3.1f]',2*pi*sample_dz*gammaHz*Ge(3),pi/1,Ge(3)));
end;

% Verify Acquisition bw
if (abs(1/dta) - (fb) > 1e-8)
	disp('Warning: Acquisition sampling rate is too low to span acq bandwidth');
	disp('Warning: Acquisition sampling rate is too low to span acq bandwidth');
% 	uiwait(msgbox('Error: Acquisition sampling rate is too low to span acq bandwidth'));
% 	error(0);
end;


