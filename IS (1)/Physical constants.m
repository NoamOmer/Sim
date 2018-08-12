
% - - - - - - - - - -
% Physical constants
% - - - - - - - - - -

global e;
global c;
global m_e;
global h_bar;
global sigma_x;
global sigma_y;
global sigma_z;
global S_x;
global S_y;
global S_z;
global I
global meu;

e     = 1.6022E-19;                % [coulomb]
c     = 3E+8;                      % [m/sec]
m_e   = 9.1095E-31;                % [kg]
h_bar = (6.6262E-34)/(2*pi);       % [J*sec] = [kg*m^2/sec]

sigma_x = [0, 1; 1 0];             % [none] Pauli matrices, dimensionless
sigma_y = [0,-i; i 0];
sigma_z = [1, 0; 0 1];
S_x     = (1/2)*h_bar*sigma_x;     % [angular momentum] Spin (operators in QM)  { =[J*sec]=[h_bar]=[kg*m^2/sec] }
S_y     = (1/2)*h_bar*sigma_y;
S_z     = (1/2)*h_bar*sigma_z;
I       = [S_x,S_y,S_z];

meu = gamma_T*I;                     % [amper*m^2] magnetic moment  { [amper*sec/kg]X[kg*m^2/sec] = amper*m^2 }

