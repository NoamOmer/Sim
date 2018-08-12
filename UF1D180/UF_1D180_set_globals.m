
set_globals;
global sim_name;         sim_name         = 'UF 1D180 Simulation';

global SIMULATION;       SIMULATION        = 1;                  % A flag: mark this run for simulation only
global Inhomo_Flag;      Inhomo_Flag       = 1;                  % A flag: enable/disable B0 inhomogeneity
global Inhomo_Fix_Flag;  Inhomo_Fix_Flag   = 0;                  % A flag: enable/disabe B0 inhomogeneity fix
global Map_Mult_Factor;  Map_Mult_Factor   = 0;                  % A flag: enable/disabe zero inhomogeneity test
global RFnGRAD_Win_Flag; RFnGRAD_Win_Flag  = 0;                  % A flag: enable/disabe RF & GRAD windowing
global RevExcAcq;        RevExcAcq         = 0;                  % A flag: Reversed (1) or equal (0) excitation
                                                                 %         & acquisition directions

global M0;               M0     = [0,0,1];                       % [amper/m]   Magnetization
global M_init;
global B0_lab;           B0_lab = [0,0,7];                       % [Tesla]     Lab magnetic field

global omega_0;          omega_0   = gamma_T*B0_lab(3);          % [rad]       Larmor frequency
global omega_rot;        omega_rot = omega_0;                    % [rad]       Rotating frame
global omega_cs;         omega_cs  = 0;                          % [rad]       Chemical shift
global tof;
global RF_ESP;           RF_ESP           = 24;
global RF_ESP_SLR;       RF_ESP_SLR       = 16;
global Ga_ESP;           Ga_ESP           = 140;

% Inhomogeneity  %
global P_DNu0;
global dOmega_0;                                                 % [Hz] % dOmega_0 = Eta0+Eta1*z_axis+Eta2*(z_axis.^2)+Eta3*(z_axis.^3);
global dB0z;                                                     % [T]

% pi/2 excitation
global B1MaxampLabHz;    B1MaxampLabHz = 16E+3;                  % [Hz]        RF field lab maximal amplitude
global dt;               dt = (31.25/16)*1e-6;                   % [sec]       Temporal resolution (during excitation)

% Sample axis
global sample_Zi;        sample_Zi  = -2.25;                     % [cm]      Sample z-axis
global sample_Zf;        sample_Zf  = +2.25;
global sample_Lz;        sample_Lz  = abs(sample_Zi-sample_Zf);
global sample_z_axis;

% 'Z' Excitation axis
global Zie;              Zie  = -1.5;
global Zfe;              Zfe  = +1.5;
global Zia;
global Zfa;
global Lz;                                                       % [cm]        Sample length
global dz;                                                       % [cm]        Spatial resolution
global z_axis;                                                   % [cm]
global Ze_INTERP_FACTOR; Ze_INTERP_FACTOR = 3.00;                %      Interpolate z-axis for higher spatial res.

% Purge %
global Tpr;              Tpr = 1000E-6;                          % [sec]       Purge duration
global Gpr;                                                      % [G/cm]      Purge gradient
global tpr;              tpr = 0:dt:Tpr;

% Excitation %
global Tp;               Tp  = 2E-3;                             % [sec]       Duration of Chirp pulse
global dte;                                                      % [sec]       Temporal resolution
global te;                                                       % [sec]       Time axis
global te_of_z;                                                  % [sec]       Time axis vs. spatial location 'z'
global O_of_te;                                                  % [Hz]        O(kHz)
global Oi;                                                       % [Hz]        O(kHz)
global R_of_te;                                                  % [1/sec^2]   O(kHz/ms)
global DOkHz;                                                    % [kHz]       Excitation spectral range
% B_RF amplitude in rot. frame: on the 300MHz the maximal value is ~8kHz
global B1ampLabHz;                                               % [Hz]        Non time-modulated RF amp O(kHz)
global B1ampCoeff;                                               % [none]      B1 amp coefficient
global Ge;               Ge  = +5;                               % [G/cm]      Magnetic field gradient
% Ge, te & O are set by inhomogeneity correction algorithm
B1ampCoeff = 0.52*3.5; % Determined manually (just run the simulation and check at what value the excitation is exactly pi/2)

%  Acquisition  %
global dZa;              %dZa = 0.5e-1;                          % [cm]   Acquisition resolution (pixel size)
global hi_res_za;                                                % [cm]   High res. acq. z-axis (correlating to 'ta')

global Ta;               Ta  = 10E-3;                            % [sec]  Duration of Chirp pulse
global ta;                                                       % [sec]  Time axis vs. hi_res_za    (Zia-->Zfa)
global ta_of_z;                                                  % [sec]  Time axis vs. basic z_axis (Zie-->Zfe)
global dta;                                                      % [sec]  Temporal resolution
global Ga_of_ta;                                                 % [G/cm] Magnetic field gradient

% global theoretical_phi_pre_acq1;
% global theoretical_phi_pre_acq2;
% global theoretical_phi_pre_acq3;

