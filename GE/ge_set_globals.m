
set_globals;

% Flags
global Inhomo_Flag;      Inhomo_Flag = 0;                                     % A flag: enable/disable B0 inhomogeneity

%  ===================
%   Manual parameters
%  ===================
global sample_Nz;        sample_Nz        = 2e4;
global required_res;     required_res     = 0.10;                             % [cm]
global Ta;               Ta               = 3e-3;                             % [sec]      Acquisition duration
Lhalf = 4/2;                                                               % [cm]

% -------------
%  Image axis
% -------------
global Zie;              Zie              = -Lhalf;                           % [cm] Image z-axis beginning
global Zfe;              Zfe              = +Lhalf;                           % [cm] Image z-axis end

% -------------
% Excitation
% -------------
global flip_a;           flip_a           = 15;                               % [deg]      Flip angle
global Gss;                                                                   % [G/cm]     Excitation gradient
global dte;              dte              = 5e-6;                             % [sec]      Temporal resolution (during excitation)
global tof;              tof              = 0;                                % [Hz]       Reference frequency

% -------------
% Acquisition
% -------------
global Dk;               Dk               = 1/required_res;                   % [cm]
global Ga;               Ga               = Dk / (gammaHz*Ta);                % [G/cm]     Acquisition gradient
global fb;               fb               = gammaHz*Ga*Lhalf*2;               % [Hz]       filter bandwidth
global dta;              dta              = 1.0*1/fb;                         % [sec]      Acquisition dwell time
global N_acq_pts;        N_acq_pts        = ceil(Ta/dta);
global dza;              dza              = Lhalf*2/N_acq_pts;
global za_axis;          za_axis          = linspace(-Lhalf+dza/2,+Lhalf-dza/2,N_acq_pts);
global ta;               ta               = linspace(0,Ta,N_acq_pts);

% -------------
% Purge
% -------------
global Tpr;              Tpr              = Ta/2;                             % [sec]      Purge duration
global Gpr;              Gpr              = -Ga;                              % [G/cm]     Magnetic field gradient

% -------------
% TR
% -------------
global TR;               TR               = 10e-3;                            % [sec]
global nTRs;             nTRs             = 3;                              % [none]


%  -------------
%   Sample axis
%  -------------
Shalf = Lhalf;
global sample_Zi;        sample_Zi        = +Shalf*sign(Zie);                 % [cm]      Sample z-axis
global sample_Zf;        sample_Zf        = +Shalf*sign(Zfe);
global sample_Lz;        sample_Lz        = Shalf*2;
global sample_dz;        sample_dz        = sample_Lz / sample_Nz;
global sample_z_axis;    sample_z_axis    = linspace(sample_Zi, sample_Zf, sample_Nz);

% -------------
%  Misc
% -------------
% Lab magnetic field
% global B1MaxampLabHz;    B1MaxampLabHz    = 16E+3;                          % [Hz]       RF field lab maximal amplitude
global B1MaxampLabHz;    B1MaxampLabHz    = 2*22e-6*gamma_Hz_T;               % [Hz]       Siemens
 
% Initial magnetization
global M0;               M0               = [0,0,1];                          % [amper/m] Magnetization

%  Inhomogeneity  %
global P_DNu0;
global dOmega_0;                                                              % [Hz]
global dB0z;                                                                  % [T]

% Frequency axis
global h2o_nu;           h2o_nu  = 0;                                         % [Hz]  Water frequency
global max_nu;           max_nu  = fb/2;                                      % [Hz]  Maximal sampled frequency
global dnu;              dnu     = 1/Ta;                                      % [Hz]  Frequency axis resolution
global nu_axis;          nu_axis = (-max_nu:dnu:max_nu-dnu)+h2o_nu;           % [Hz]  Frequency axis

% Verify simulation spatial axis
if (sample_dz*gammaHz*Ga > pi/3)  % pi/3 instead of pi or 2*pi just to be on the safe side.
	error(sprintf('sample spatial resolution is too low (sample_dz*gammaHz*Ga=%3.2f)',sample_dz*gammaHz*Ga));
end;

% Verify Acquisition bw
if (round(1/dta) < round(abs(fb)))
	uiwait(msgbox('Error: Acquisition sampling rate is too low to span acq bandwidth'));
	error(0);
end;
