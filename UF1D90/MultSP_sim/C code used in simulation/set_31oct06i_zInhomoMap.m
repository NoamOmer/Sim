% at  =     1.50[ms] 	NTau    =    71.00     	Gz_p    =    -9.00[G/cm]
% sw  =   100.00[kHz]	tpwr    =    56.00[dB] 	Gz_a    =     9.00[G/cm]
% np  =   300.00     	gain    =    40.00[dB]	Tau     =     0.00[ms]  
% fb  =    40.00[kHz]	presig  =    l      	dTau    =     4.00[ms]  
% dt  =       10[us] 	t_p     =     0.75[ms]	tof_A   =   484.00[Hz]  
% nt  =        1     	t_delay =    10.00[us]	d1      =    15.00[sec] 
% ss  =        3     	dNu     =   3.5714[Hz] 	dz      =   0.1740[mm]	
% ct  =        1     	DeltaNu = 250.0000[Hz] 	Delta_z =  26.1008[mm]	
% Gname = gwurst

gamma_T = 4.257;           % [kHz]
at       = 1.5;             % [ms]
Ta       = at;              % [ms]
np       = 300;
Nk       = np/2;
dt       = 10E-3;           % [ms]
NTau     = 71;
dTau     = 4;               % [ms]
Gz_a     = 9;               % [G/cm]

DEFAULT_CENTER_POINT = round(Nk/2);