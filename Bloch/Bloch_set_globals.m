
set_globals;

global M0;                M0 = (1/2)*[0,0,1];              % [amper/m]   Magnetization

global B0_lab;            B0_lab = [0,0,1];                % [Tesla]     Lab magnetic field
global B1MaxampLabHz;     B1MaxampLabHz = 0.1E+3;          % [Hz]        RF field lab amplitude
global omega_0;           omega_0   = gamma_T*B0_lab(3);   % [rad]       Larmor frequency
global omega_rot;         omega_rot = omega_0;             % [rad]       Rotating frame

global dt;                dt = 0.02E-3;                    % [sec]       Temporal resolution
global T;                 T  = 10E-3;                      % [sec]
global t;                 t  = 0:dt:T+dt;                  % [sec]       Time axis]

