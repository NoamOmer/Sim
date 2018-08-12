
set_globals;
global sim_name;         sim_name         = 'UF 1D Simulation';
global Seq1D;            Seq1D            = 1; 
global Seq2D;            Seq2D            = 0;

global SIMULATION;       SIMULATION       = 1;                   % A flag: mark this run for simulation only
global Inhomo_Flag;      Inhomo_Flag      = 0;                   % A flag: enable/disable B0 inhomogeneity
global DesignPulse_Flag; DesignPulse_Flag = 0;                   % A flag: enable/disabeB0 inhomogeneity fix using designed pulses
global Inhomo_Fix_Flag;  Inhomo_Fix_Flag  = 0;                   % A flag: enable/disabeB0 inhomogeneity fix
global RFnGRAD_Win_Flag; RFnGRAD_Win_Flag = 0;                   % A flag: enable/disabe designed Chirp windowing
global SLR_flag;         SLR_flag         = 1;                   % A flag: enable/disabe SLR pulse
global SE_flag;          SE_flag          = 0;                   % A flag: enable/disabe SE PI pulse

global M0;               M0               = [0,0,1];             % [amper/m] Magnetization
global M_init;
global B0_lab;           B0_lab           = [0,0,7];             % [Tesla]   Lab magnetic field

global omega_0;          omega_0          = gamma_T*B0_lab(3);   % [rad]     Larmor frequency
global omega_rot;        omega_rot        = omega_0;             % [rad]     Rotating frame
global omega_cs;         omega_cs         = 0;                   % [rad]     Chemical shift
global RF_ESP;           RF_ESP           = 16;
global RF_ESP_SLR;       RF_ESP_SLR       = 16;
global Ga_ESP;           Ga_ESP           = 140;

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

%  ---------------------
%   'Z' EXCITATION axis
%  ---------------------
global Zie;              Zie              = -1.75;               % [cm] Excitation z-axis beginning
global Zfe;              Zfe              = +1.75;               % [cm] Excitation z-axis end
global Zia;
global Zfa;
global Lz;               Lz               = abs(Zfe-Zie);        % [cm] Excitation length
global dz;                                                       % [cm] Spatial resolution
global z_axis;                                                   % [cm] EXCITATION (!) z-axis 
global Ze_INTERP_FACTOR; Ze_INTERP_FACTOR = 3.0;                 %      Interpolate z-axis for higher spatial res.
global alpha1;
global alpha2;
                                                                 %      Zie and Zfe will fit into it
%  -------------
%   Sample axis
%  -------------
global sample_Zi;        sample_Zi        = +2.00*sign(Zie);     % [cm]      Sample z-axis
global sample_Zf;        sample_Zf        = +2.00*sign(Zfe);
global sample_Lz;        sample_Lz        = abs(sample_Zi-sample_Zf);
global sample_z_axis;

%  ---------------
%   Inhomogeneity
%  ---------------
global Map_Mult_Factor;  Map_Mult_Factor  = 0;                   % A flag: enable/disabe zero inhomogeneity test
global gPolyOrder;       gPolyOrder       = 4;
global P_DNu0;
global dOmega_0;                                                 % [Hz] dOmega_0 = Eta0+Eta1*z_axis+Eta2*(z_axis.^2)+Eta3*(z_axis.^3);
global dB0z;                                                     % [T]

%  ------------
%   Excitation
%  ------------
global Tp;               Tp               = 3E-3/Nsp;            % [sec]     Duration of Chirp pulse
global rfwdth;           rfwdth           = Tp + 2E-3;           % [sec]
global Ge;               Ge               = [0,0,10];            % [G/cm]    Magnetic field gradient
global dte;                                                      % [sec]     Temporal resolution
global te;                                                       % [sec]     Time axis
global te_of_z;                                                  % [sec]     Time axis vs. spatial location 'z'
global O_of_te;                                                  % [Hz]      O(kHz)
global Oi;                                                       % [Hz]      O(kHz)
global R_of_te;                                                  % [1/sec^2] O(kHz/ms)
global DOkHz;                                                    % [kHz]     Excitation spectral range

global OmegaE;
global ze_of_OmegaE;
global Phi_e;                                                    % [none]    Phase at end of excitation

% B_RF amplitude in rot. frame: on the 300MHz the maximal value is 8kHz
global B1MaxampLabHz;  B1MaxampLabHz = 16E+3;                    % [Hz]      RF field lab maximal amplitude
global B1ampLabHz;                                               % [Hz]      Non time-modulated RF amp O(kHz)
global B1ampCoeff;                                               % [none]    B1 amp coefficient
if (DesignPulse_Flag)
	% Ge, te & O are set by inhomogeneity correction algorithm
	B1ampCoeff = 0.53; % Determined manually (just run the simulation and check at what value the excitation is exactly pi/2)
else
	dte = 5E-6;
	te  = 0:dte:(rfwdth-dte);
	if (~isempty(Lz))
		Oi  = gammaHz*Ge(3)*Zie;
		Of  = gammaHz*Ge(3)*Zfe;
		dO  = (Of - Oi)/Nsp;
		DOkHz = 1E-3*dO;
		R_of_te = dO / Tp;
		O_of_te = linspace(Oi,Of,length(te)); %Oi + R_of_te*te;
		DO_Hz = abs(O_of_te(end)-O_of_te(1));
	end;
	% 0.26 constant is given in JMR 172 (2005) 179–190.
	% We use 0.52 for the lab frame and then divide it by two for the rotating frame.
	B1ampCoeff = 0.52;
	B1ampLabHz = B1ampCoeff*sqrt(abs(R_of_te));
end;

%  -------------
%   Acquisition
%  -------------
% 'Z' acquisition axis
global dZa;                                                      % [cm]   Effective acquisition resolution
global za_axis;                                                  % [cm]   Acquisition z-axis
global hi_res_za;                                                % [cm]   High res. acq. z-axis (matching 'ta')

global Ta;               Ta               = Tp; %4E-3/Nsp;            % [sec]  Acquisition duration
global ta;                                                       % [sec]  Time axis vs. hi_res_za    (Zia-->Zfa)
global ta_of_z;                                                  % [sec]  Time axis vs. basic z_axis (Zie-->Zfe)
global dta;                                                      % [sec]  Temporal resolution
global GRADdt;           GRADdt           = 4e-6;                % [sec]  Acquisition temporal resolution
global Ga;                                                       % [G/cm] Acquisition gradient
global Ga_of_ta;                                                 % [G/cm] Time-dependent acquisition gradient
if (DesignPulse_Flag)
	% ta, dta & Ga are set by inhomogeneity correction algorithm
else
	dta = 50E-6;
	ta  = 0:dta:(Ta-dta);
	Ga  = (Ge*Tp/Ta) * sign(SE_flag-0.5);    % reverse Ga in case a SE was performed
	dZa = sqrt((Lz/Nsp)/(gammaHz*Ge(3)*Tp));
end;

