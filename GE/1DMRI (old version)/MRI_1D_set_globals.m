
set_globals;

global M0;                M0 = [0,0,1];                       % [amper/m]   Magnetization

global Lz;                Lz = 2.0E-0;                        % [m]         Sample length
global dz;                dz = Lz/1202;                       % [m]         Spatial resolution
global z_axis;            z_axis = (-Lz/2):dz:(Lz/2-dz);

global B0_lab;            B0_lab = [0,0,7];                   % [Tesla]     Lab magnetic field
global omega_0;           omega_0   = gamma_T*B0_lab(3);      % [rad]       Larmor frequency
global omega_rot;         omega_rot = omega_0;                % [rad]       Rotating frame
global omega_cs;          omega_cs = 0;                       % [rad]       Chemical Shift

global B1MaxampLabHz;     B1MaxampLabHz = 16E+3;              % [Hz]        RF field lab maximal amplitude
global GRAD;              GRAD = [0 0 0.03];                  % [T/m]       Magnetic field gradient

% dt < (2*pi)/(gamma_T*GRAD(3)*Lz)
global dt;                dt = 0.02E-3;                       % [sec]       Temporal resolution
global dk;                dk = gamma_T*GRAD(3)*dt;            % [rad/m]     k-space resolution

%global T_grad;           T_grad = 0.3;                       % [sec]       Gradient application time
global T_grad;            T_grad = round(((2*pi)/(gamma_T*GRAD(3)*dz))/dt)*dt;

