
set_globals;
global root_dir;       root_dir    = 'D:\PhD\06 Experiments_1D\1D UF Self Refocusing MRI\SeqParams\';

% Flags
global SLR_flag;       SLR_flag    = 1;
global se_flag;        se_flag     = 0;
global orient_flag;    orient_flag = 1;
global ovs_flag;       ovs_flag    = 0;

% General definitions
global tep;            tep         = 35;                                        % [us]
global tof;            tof         = -1052;                                     % [Hz]
global gain;           gain        = 36;                                        % [dB]
global tDel;           tDel        = 4E-6;                                      % [sec]
global sw;             sw          = 250e+3;                                    % [Hz]

% Spatial axis
global Zie;            Zie         = -3.75;                                     % [cm]
global Zfe;            Zfe         = +3.75;                                     % [cm]
global Lz;             Lz          = abs(Zie-Zfe);                              % [cm]

% Excitation
global Nsp;            Nsp         = 2;

global tpwr90;         tpwr90      = 60;                                        % [dB]
global pw90;           pw90        = 23E-6;                                     % [sec]

global Ge;             Ge          = 7.8;                                       % [G/cm]
global Tp;             Tp          = 6e-3/Nsp;                                  % [sec]
global rfwdth;         rfwdth      = Tp + 2E-3;                                 % [sec]
global dte;            dte         = 3e-6;                                      % [sec]
global rf_R;           rf_R        = gammaHz*Ge*(Lz/Nsp)/Tp;                    % [Hz/sec]
global nRF;            nRF         = round(rfwdth/dte);                         % [none]
global Gename;         Gename      = 'nbe_T5000G10_rt20ft20c4_96L_tCF1p0048';

global z_axis;         z_axis      = linspace(Zie,Zfe,nRF);                     % [cm]
global dz;             dz          = Lz/nRF;                                    % [cm]

global RF_ESP_SLR;     RF_ESP_SLR  = 32;

global RFPhaseSign;
if (SLR_flag)
	RFPhaseSign = +1;
else
	RFPhaseSign = -1;
end;

% Acquisition
global Ta;             Ta          = Tp;                                      % [sec]
global dta;            dta         = 1/sw;                                    % [sec]
global Ga;             Ga          = -Ge;                                     % [G/cm]
global Ganame;         Ganame      = 'nbe_T5000G10_rt20ft20c4_96L_tCF1p0048';

