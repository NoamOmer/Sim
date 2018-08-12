
sample_Nz        = 8000;

Lhalf            = 2.0/2;                            % [cm]
Zie              = -Lhalf;                           % [cm] Image z-axis beginning
Zfe              = +Lhalf;                           % [cm] Image z-axis end

Shalf = Lhalf;
sample_Zi        = +Shalf*sign(Zie);                 % [cm]      Sample z-axis
sample_Zf        = +Shalf*sign(Zfe);
sample_Lz        = Shalf*2;
sample_dz        = sample_Lz / sample_Nz;
sample_z_axis    = linspace(sample_Zi, sample_Zf, sample_Nz);

sample_z_axis_generic    = linspace(-1, 1, 1e6);
sample1          = exp(-8*((sample_z_axis_generic).^16));

sample1 = sample1(round(linspace(1,length(sample1),sample_Nz)));
fp(sample1);

