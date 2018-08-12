% - - - - - - - - - - - - - - -
% Global path Parameters
% - - - - - - - - - - - - - - -
global	OS_COPY;
% if (ismac)
	OS_COPY  = 'cp';


	global extRFpulseDir;           extRFpulseDir         = '/Users/noambe/Dropbox/E/Siemens/04 extRF Pulse lib/';
	global extRFpulseDir_Bruker;    extRFpulseDir_Bruker  = '/Users/noambe/OneDrive - Tel-Aviv University/MRI_data/Bruker_7T/Bruker_sources/';
	global EMC_root;                EMC_root               = nbePath('/Users/noambe/Dropbox/EMC_T2_FIT/');
% else
% 	OS_COPY  = 'copy';
% 	DriveDir = 'E:\';
% 	global extRFpulseDir;           extRFpulseDir          = nbePath(['C:\Users\beneln01\Dropbox\E\01 Post\06 Projects\Siemens\04 extRF Pulse lib\']);
% 	global extRFpulseDir_Bruker;    extRFpulseDir_Bruker  = '/Users/noambe/OneDrive - Tel-Aviv University/MRI_data/Bruker_7T/Bruker_sources/';
% 	global SimRoot;                 SimRoot                = nbePath(['C:\Users\beneln01\Dropbox\1\Sim\']);
% 	global EMC_root;                EMC_root               = nbePath(['C:\users\beneln01\Dropbox\EMC_T2_FIT\']);
% end;

% - - - - - - - - - - - - - - -
% Global Simulation Parameters
% - - - - - - - - - - - - - - -
% global interactive_mode; interactive_mode = 1;

% warning off;
warning('off', 'all');

global DEBUG_FLAG;
% 0, 1, 2, 3, 4 or 5
% if isempty(DEBUG_FLAG)
	DEBUG_FLAG = 4;
% end;
global SIMULATION;  SIMULATION = 1;                   % A flag: mark this run for simulation only

global gammaHz;     gammaHz    = 2.675222005*1e8/(1e4*2*pi);   % [Hz/G]
global gammakHz;    gammakHz   = gammaHz/1000;                 % [kHz/G]
global gammaC13Hz;  gammaC13Hz = gammaHz*79.4/301;             % [Hz/G]
global gamma_T;     gamma_T    = gammaHz*2*pi*1e+4;            % [rad/(sec*Tesla)] = [rad*amper*sec/kg]  ('H' protons)
global gamma_Hz_T;  gamma_Hz_T = gammaHz*1e+4;                 % [Hz/Tesla]

global h_bar;       h_bar      = 6.62607004e-34;               % [Joule*sec]
global kB;          kB         = 1.38064852e-23;               % [Joule/Kelvin]
global Avogadro_C;  Avogadro_C = 6.02214086e+23;               % [mol-1]

% Rotating frame parameters
global B_rot;       B_rot      = [0,0,0];
global omega_0;     omega_0    = 0;
global omega_rot;   omega_rot  = 0;
global omega_CS;    omega_CS   = 0;

global GoldenAngle; GoldenAngle= 111.246117975;       % [deg] Golden angle

global Grs_lmt;                                       % [G/(cm*sec)] Gradient rise time limit
if SIMULATION 
	Grs_lmt = 5*1e6;
else
	Grs_lmt = 0.66*1e6; % 660*1e3
end;


global B1MaxampLabHz;    B1MaxampLabHz    = 2*22e-6*gamma_Hz_T;               % [Hz]       Siemens


% Chirp pulse calibration
global B1Chirp90Coeff;   B1Chirp90Coeff  = 0.265;
global B1Chirp180Coeff;  B1Chirp180Coeff = 0.265*3.3;


global RH_flag;     RH_flag    = 1;                   % A flag: Right- or Left- Hand-Rule rotation

global Relax_Flag;  Relax_Flag = 1;                   % Activates / deactivates relaxation
% global T2;          T2         =  60E-3;             % [sec]
% global T1;          T1         = 100E-3;             % [sec]

global Rect_P;      Rect_P     = 1;
global Sinc_P;      Sinc_P     = 2;

% plot styles
global styles;      styles = {'k-','b-','r-','g-','m-','k.-','b.-','r.-','g.-','m.-','k--','b--','r--','c--','m--','k.-','b-','r.-','g-','m.-','k.-','b--','r.-','g-','m.-','k-','b.-','r-','c.-','m-'};

% Run magnetization evolution simulations using C code
global C_CODE;      C_CODE     = 1;

% % 2D pulses - Minimal number of gradient plateau points
% global MIN_2D_PLATEAU_PTS;  MIN_2D_PLATEAU_PTS  = 1;
% global K_SPACE_TRAJ_ZIGZAG; K_SPACE_TRAJ_ZIGZAG = 1;
% global K_SPACE_TRAJ_CARTS;  K_SPACE_TRAJ_CARTS  = 2;
