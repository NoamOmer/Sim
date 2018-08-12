% - - - - - - - - - - - - - - -
% Global Simulation Parameters
% - - - - - - - - - - - - - - -

global gamma_T;       gamma_T   = 2.675E+8;                   % [rad/(sec*Tesla)] = [rad*amper*sec/kg]  ('H' protons)
global M0;            M0 = (1/2)*[0,0,1];                     % [amper/m]   Magnetization

global omega_0;       omega_0   = gamma_T*B0_lab(3);          % [rad]       Larmor frequency
global omega_rot;     omega_rot = omega_0;                    % [rad]       Rotating frame
global omega_cs;      omega_cs = 0;                           % [<<>>]      Chemical shift

global B0_lab;    B0_lab = [0,0,1];                           % [Tesla]     Lab magnetic field
global B0_rot;    B0_rot = B0_lab-[0,0,(omega_rot/gamma_T)];  % [Tesla]     Rotating frame magnetic field

global Lz;            Lz = 6.0E-2;                            % [m]         Sample length
global dz;            dz = 0.1E-3;                            % [m]         Spatial resolution
global z_axis;        z_axis = (-Lz/2):dz:(Lz/2-dz);

global Tp;            Tp = 10E-3;                             % [sec]
global dt;            dt = 0.02E-3;                           % [sec]       Temporal resolution
global t;             t  = 0:dt:T;

%   C h i r p   %
global B1MaxampLabHz; B1MaxampLabHz = 0.1E+3;                 % [Hz]        RF field lab maximal amplitude
global Oi;            Oi = 1E+3;                              % [Hz]        O(kHz)
global R;             R   = 1E+6;                             % [Hz/sec]    O(kHz/ms)

%   G r a d i e n t s   %
global Ge;            Ge = [0 0 0.001];                       % [T/m]       Magnetic field gradient
global G_acq;         G_acq = [0 0 0.001];                    % [T/m]       Magnetic field gradient






global gamma_T;         gamma_T   = 2.675E+8;                  % [rad/(sec*Tesla)] = [rad*amper*sec/kg]  ('H' protons)

global B0_lab;          B0_lab = [0,0,1];                      % [Tesla]     Lab magnetic field
global B1MaxampLabHz;   B1MaxampLabHz = 0.1E+3;                % [Hz]        RF field lab maximal amplitude
global GRAD;            GRAD = [0 0 0.001];                    % [T/m]       Magnetic field gradient

% dt < (2*pi)/(gamma_T*GRAD(3)*Lz)
global dt;              dt = 0.02E-3;                          % [sec]       Temporal resolution
global dk;              dk = gamma_T*GRAD(3)*dt;               % [rad/m]     k-space resolution
%global T_grad;          T_grad = 0.3;                          % [sec]       Gradient application time
global T_grad;          T_grad = round(((2*pi)/(gamma_T*GRAD(3)*dz))/dt)*dt;

global omega_0;         omega_0   = gamma_T*B0_lab(3);         % [rad]       Larmor frequency
global omega_rot;       omega_rot = omega_0;                   % [rad]       Rotating frame

global Relax_Flag;      Relax_Flag = 0;                        % Activates / deactivates relaxation
global T2;              T2 = 50E-3;                            % [sec]
global T1;              T1 = 500E-3;                           % [sec]

global Rect_P;          Rect_P = 1;
global Sinc_P;          Sinc_P = 2;

