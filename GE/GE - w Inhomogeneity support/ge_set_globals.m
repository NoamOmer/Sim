
set_globals;

% Flags
global Inhomo_Flag;      Inhomo_Flag = 0;                          % A flag: enable/disable B0 inhomogeneity

% Sample
global M0;               M0        = [0,0,1];                      % [amper/m]  Magnetization

% Lab magnetic field
global B1MaxampLabHz;    B1MaxampLabHz = 16E+3;                    % [Hz]       RF field lab maximal amplitude
global dt;               dt = (31.25/16)*1e-6;                     % [sec]      Temporal resolution (during excitation)

global B0_lab;           B0_lab    = [0,0,7];                      % [Tesla]    Lab magnetic field
global B_rot;            B_rot     = [0,0,0];                      % [Tesla]    Rotating frame magnetic field
global omega_0;          omega_0   = gamma_T*B0_lab(3);            % [rad]      Larmor frequency
global omega_rot;        omega_rot = omega_0;                      % [rad]      Rotating frame
global omega_cs;         omega_cs  = 0;                            % [rad]      Chemical shift
global tof;              tof       = -558.3;                       % [Hz]       Reference frequency

%  Inhomogeneity  %
global P_DNu0;
global dOmega_0;                                                   % [Hz]
global dB0z;                                                       % [T]

% Excitation - Range should be set according to the maximal FOV in the Inhomogeneity map
global flip_a;           flip_a = 90;                              % [deg]      Flip angle
global Zie;              Zie    = -1.83;                           % [cm]       Sample start point
global Zfe;              Zfe    = +1.93;                           % [cm]       Sample end point
global Lz;               Lz     = abs(Zfe-Zie);                    % [cm]       FOV
global dz;               dz     = 1E-3;                            % [cm]       Spatial resolution
global z_axis;                                                     % [cm]       Spatial axis
global Z_EXTRAP_FACTOR ; Z_EXTRAP_FACTOR = 1.5;                    % [none]     Spatial axis extrapolation factor

%  Acquisition  %
global Ta;               Ta    = 100e-3;                           % [sec]      Acquisition duration
global sw;               sw    = 100e+3;                           % [Hz]       Spectral width
global dta;              dta   = 1/sw;                             % [sec]      Acquisition temporal resolution
global n_pts;            n_pts = 2*Ta/dta;  % NOTE: # COMPLEX PTS  % [none]     Number of COMPLEX acquisition points
global ta;               ta    = 0:dta:Ta-dta;                     % [sec]      Acquisition time axis
global Ga;               Ga    = 1.00;                             % [G/cm]     Magnetic field gradient

%  Purge  %
global Tpr;              Tpr   = Ta/2;                             % [sec]      Purge duration
global Gpr;              Gpr   = -Ga;                              % [G/cm]     Magnetic field gradient

% Frequency axis
global h2o_nu;           h2o_nu  = 950;                            % [Hz]  Water frequency
global max_nu;           max_nu  = sw/2;                           % [Hz]  Maximal sampled frequency
global dnu;              dnu     = 1/Ta;                           % [Hz]  Frequency axis resolution
global nu_axis;          nu_axis = (-max_nu:dnu:max_nu-dnu)+h2o_nu;% [Hz]  Frequency axis

