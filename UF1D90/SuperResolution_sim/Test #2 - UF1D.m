% 1D Spatially-Encoded MRI Simulation
clcl;

% -------------------------------------------------------------------------
% Complie C Code
% -------------------------------------------------------------------------
cur_dir = pwd;
cd('D:\PhD\Matlab\Simulations\Ver_4\C code');
mex 'D:\PhD\MATLAB\Simulations\Ver_4\C code\evolve_M_n_acq1Dgrad_t_CPP.c'
cd(cur_dir);

% -------------------------------------------------------------------------
% Set context
% -------------------------------------------------------------------------
context = 'UF_1D_set_globals';
set_context;
declare_start(sim_name);  tic
disp(sprintf('Expected dZa = %3.3f [mm]',dZa*10));

% -------------------------------------------------------------------------
% Multi-SP parameter optimization
% -------------------------------------------------------------------------
if (Nsp > 1)
	nn = 0:3;
	OPnp = round((gammaHz * Ge(3) * Tp * Lz) ./ (2*nn + 1));
	declare_stage('Multi-SP paraemters optimization');
	disp(sprintf('Optimum    np  :  %3.0f  %3.0f  %3.0f  %3.0f',OPnp(1),OPnp(2),OPnp(3),OPnp(4)));
	disp(sprintf('Simulation np  :  %3.0f                     ',Ta/dta                         ));
	
	nMz = Ta/dta;
	dx  = Lz / (2*nMz);
	Dk  = 1/dx;
	k0_tag = gammaHz*Ge(3)*Tp;
	k0_tag = mod(k0_tag,Dk);
	disp(sprintf('\n  Dk = %3.1f [%3.1f ... %3.1f]\n  k0 = %3.1f\n',Dk,-Dk/2,+Dk/2,k0_tag));
end;

% -------------------------------------------------------------------------
% Retrieve the B0 inhomogeneity
% -------------------------------------------------------------------------
[z_axis,P_DNu0,fh] = UD1D_retrieve_inhomogeneity(context);

% -------------------------------------------------------------------------
% Set the spatial axis
% -------------------------------------------------------------------------
% Interpolate the internal z_axis so that the RF dwell time in the te axis
% calculation will match the global dwell-time value
n_RF_pts  = 401; %1 + Tp/dte;
z_axis = linspace(Zie,Zfe,n_RF_pts*Ze_INTERP_FACTOR);
dz = z_axis(2) - z_axis(1);
% Lz = abs(z_axis(end) - z_axis(1));
disp(sprintf('UF1D180:\n dz=%d;\n Lz=%d',dz,Lz));
if (fh)
	plot(z_axis,polyval(P_DNu0,z_axis),'m.');
	legend({'Initial','After Inhomo scaling','After z-interp and range setting'},'Location','Best'); grid;
end;

% -------------------------------------------------------------------------
% Create a shaped sample
% -------------------------------------------------------------------------
[sample_z_axis,M_init] = UD1D_create_sample(context);

% Now that the z_axis is arranged -- define the inhomogeneity field
B_lab = B0_lab + [0, 0, (omega_cs /gamma_T)];
B_rot = B_lab  - [0, 0, (omega_rot/gamma_T)];
dOmega_0 = polyval(P_DNu0,sample_z_axis);      % [Hz]
dB0z = (dOmega_0/gammaHz)*1E-4;                % [T]

% -------------------------------------------------------------------------
% Design & Apply PI/2-Chirp pulse
% -------------------------------------------------------------------------
[chirp_rot,phi_chirp_rot,alpha0,alpha1,alpha2] = UD1D_design_pulse(context);

verify_inhom_fix_1D(context);

x_chirp_rot = real(chirp_rot);
y_chirp_rot = imag(chirp_rot);
z_chirp_rot = zeros(1,length(chirp_rot));

if (DEBUG_FLAG > 1)
	figure(93); hold on;  plot(te*1e+3, x_chirp_rot, 'k.-', te*1e+3, y_chirp_rot, 'b.-');
	title('Chirp pulse in rotating frame');      xlabel('Time [ms]');  ylabel('Magnetic field'); 
	legend({'x-coordinate','y-coordinate'});     set_gca;
	figure(94); hold on;  plot(te*1e+3, abs(chirp_rot), 'b.-');
	title('Chirp amplitude in rotating frame');  xlabel('Time [ms]');  ylabel('Amp [T]'); set_gca;
	figure(92); hold on;  plot(te*1e+3, phi_chirp_rot, 'r.-');
	title('Chirp angle in rotating frame');      xlabel('Time [ms]');  ylabel('Angle [rad]'); set_gca;
end;

declare_stage('Apply PI/2-Chirp pulse');
B_eff_rot = [transpose(x_chirp_rot),transpose(y_chirp_rot),transpose(z_chirp_rot)] + (transpose(ones(1,length(chirp_rot))))*B_rot;
if (C_CODE)
	zGrad = ones(1,length(B_eff_rot(:,1)))*Ge(3);
	st = set_evolve_M_n_acq1Dgrad_t_CPP_struct(M_init,B_eff_rot,dte,dB0z,zGrad,sample_z_axis,RH_flag,...
	                                           0,Inhomo_Flag,Relax_Flag,T1,T2,context);
	[Mx,My,Mz,dummy1,dummy2] = evolve_M_n_acq1Dgrad_t_CPP(st);
	M = [transpose(Mx),transpose(My),transpose(Mz)];
else
	[M, dummy] = evolve_M_n_acq(M_init, B_eff_rot, dB0z, Ge, dte, sample_z_axis, 0, Inhomo_Flag, 0);
end;
plot_M_n_Phase(0,M_init,M,'M(z) post PI/2 Chirp pulse','Phase post PI/2 Chirp pulse',{'\phi(z)'},context);
% M2(:,1) = M_init(:,3)/sqrt(2);
% M2(:,2) = M_init(:,3)/sqrt(2);
% M2(:,3) = M_init(:,1);
% plot_M_n_Phase(0,M_init,M2,'M(z) post PI/2 Chirp pulse','Phase post PI/2 Chirp pulse',{'\phi(z)'},context);
% M = M2;

% -------------------------------------------------------------------------
% Refocus SLR PI/2-Chirp pulse
% -------------------------------------------------------------------------
if (SLR_flag)
	declare_stage('Refocus SLR pulse');
	Trefocus     = 0.5*(rfwdth-Tp);
	dt_refocus   = dte;
	nrefocus_pts = Trefocus/dt_refocus;
	x_chirp_rot = zeros(1,nrefocus_pts);
	y_chirp_rot = zeros(1,nrefocus_pts);
	z_chirp_rot = zeros(1,nrefocus_pts);
	B_eff_rot   = [transpose(x_chirp_rot),transpose(y_chirp_rot),transpose(z_chirp_rot)];
	acq_flag    = 0;

	zGrad = ones(1,length(B_eff_rot(:,1))) * (-Ge(3));
	st = set_evolve_M_n_acq1Dgrad_t_CPP_struct(M,B_eff_rot,dt_refocus,dB0z,zGrad,sample_z_axis,RH_flag,...
											   acq_flag,Inhomo_Flag,Relax_Flag,T1,T2,context);
	[Mx,My,Mz,dummy1,dummy2] = evolve_M_n_acq1Dgrad_t_CPP(st);
	M = [transpose(Mx),transpose(My),transpose(Mz)];
	plot_M_n_Phase(0,M_init,M,'M(z) post SLR refocusing','Phase post SLR refocusing',{'\phi(z)'},context);
end;
% if (0) % Gradient-Echo TEST - in order to view the phase, post excitation
% end;

% -------------------------------------------------------------------------
% Reverse acquisition direction - purge pre PI pulse
% -------------------------------------------------------------------------
rev_flag = 0;
if (rev_flag)
	declare_stage('Pre PI Pulse purge');
	Tpurge     = 500e-6;
	dt_purge   = dte;
	nPurge_pts = Tpurge/dt_purge;
	x_chirp_rot = zeros(1,nPurge_pts);
	y_chirp_rot = zeros(1,nPurge_pts);
	z_chirp_rot = zeros(1,nPurge_pts);
	B_purge_rot = [transpose(x_chirp_rot),transpose(y_chirp_rot),transpose(z_chirp_rot)];
	zGradPurge = mean(Ga_of_ta)*Ta/Tpurge;
	zGradPurge = zGradPurge/2;             % half before and half after the PI pulse
	if (SE_flag) zGradPurge = -zGradPurge; end;
	zGradPurge = zGradPurge * ones(1,length(B_purge_rot(:,1)));
	acq_flag = 0;
	st = set_evolve_M_n_acq1Dgrad_t_CPP_struct(M,B_purge_rot,dt_purge,dB0z,zGradPurge,sample_z_axis,RH_flag,...
	                                           acq_flag,Inhomo_Flag,Relax_Flag,T1,T2,context);
	[Mx,My,Mz,dummy1,dummy2] = evolve_M_n_acq1Dgrad_t_CPP(st);
	M = [transpose(Mx),transpose(My),transpose(Mz)];
	plot_M_n_Phase(0,M_init,M,'M(z) Purge 1st part','Phase Purge 1st part',{'\phi(z)'},context);

	% Reverse the acquisition gradient
	Ga_of_ta = -Ga_of_ta;
end;

% -------------------------------------------------------------------------
% Spin-Echo pulse
% -------------------------------------------------------------------------
if (SE_flag)
	declare_stage('Pre PI Pulse Crusher');
	Tcrush     = 1e-3;
	dt_crush   = dte;
	nCrush_pts = Tcrush/dt_crush;
	x_chirp_rot = zeros(1,nCrush_pts);
	y_chirp_rot = zeros(1,nCrush_pts);
	z_chirp_rot = zeros(1,nCrush_pts);
	B_crush_rot = [transpose(x_chirp_rot),transpose(y_chirp_rot),transpose(z_chirp_rot)];
	zGradCrush = 10;
	zGradCrush = zGradCrush * ones(1,length(B_crush_rot(:,1)));

	% Pre PI crusher
	st = set_evolve_M_n_acq1Dgrad_t_CPP_struct(M,B_crush_rot,dt_crush,dB0z,zGradCrush,sample_z_axis,RH_flag,...
	                                           0,Inhomo_Flag,Relax_Flag,T1,T2,context);
	[Mx,My,Mz,sig_real,sig_imag] = evolve_M_n_acq1Dgrad_t_CPP(st);
	M = [transpose(Mx),transpose(My),transpose(Mz)];

	% PI pulse
	declare_stage('PI Pulse (non slice--selective hard pulse)');
	[B_RF_rot, B_RF_t] = gen_RF_pulse(Rect_P, pi, 0, dte, omega_cs, context);
	disp(sprintf('PI Pulse duration = %d [us]',(B_RF_t(end))*1E+6));
	B_eff_rot = B_RF_rot + (transpose(ones(1,length(B_RF_t)))) * B_rot;
	acq_flag = 0;

	if (C_CODE)
		zGrad = zeros(1,length(B_RF_rot(:,1)));
		st = set_evolve_M_n_acq1Dgrad_t_CPP_struct(M,B_eff_rot,dte,dB0z,zGrad,sample_z_axis,RH_flag,...
		                                           acq_flag,Inhomo_Flag,Relax_Flag,T1,T2,context);
		[Mx,My,Mz,dummy1,dummy2] = evolve_M_n_acq1Dgrad_t_CPP(st);
		M = [transpose(Mx),transpose(My),transpose(Mz)];
	else
		[M, dummy] = evolve_M_n_acq(M,B_eff_rot,dB0z,[0,0,0],dte,sample_z_axis,acq_flag,Inhomo_Flag,DEBUG_FLAG);
	end;

	% Post PI crusher
	declare_stage('Post PI Pulse Crusher');
	st = set_evolve_M_n_acq1Dgrad_t_CPP_struct(M,B_crush_rot,dt_crush,dB0z,zGradCrush,sample_z_axis,RH_flag,...
	                                           0,Inhomo_Flag,Relax_Flag,T1,T2,context);
	[Mx,My,Mz,sig_real,sig_imag] = evolve_M_n_acq1Dgrad_t_CPP(st);
	M = [transpose(Mx),transpose(My),transpose(Mz)];

	plot_M_n_Phase(0,M_init,M,'M(z) after PI pulse','Phase after PI pulse',{'\phi(z)'},context);
end;

% -------------------------------------------------------------------------
% Reverse acquisition direction - purge post PI pulse
% -------------------------------------------------------------------------
if (rev_flag)
	declare_stage('Post PI Pulse purge');
	if (SE_flag) zGradPurge = -zGradPurge;   end;     % reverse the purge gradient sign
	acq_flag = 0;
	st = set_evolve_M_n_acq1Dgrad_t_CPP_struct(M,B_purge_rot,dt_purge,dB0z,zGradPurge,sample_z_axis,RH_flag,...
	                                           acq_flag,Inhomo_Flag,Relax_Flag,T1,T2,context);
	[Mx,My,Mz,dummy1,dummy2] = evolve_M_n_acq1Dgrad_t_CPP(st);
	M = [transpose(Mx),transpose(My),transpose(Mz)];
	plot_M_n_Phase(0,M_init,M,'M(z) Purge 2nd part','Phase Purge 2nd part',{'\phi(z)'},context);
end;

% -------------------------------------------------------------------------
% Acquire signal
% -------------------------------------------------------------------------
declare_stage('Acquiring signal');
B_eff_rot = transpose(ones(1,length(ta))) * B_rot;
acq_flag = 1;
if (DesignPulse_Flag)
	if (C_CODE)
		if (DEBUG_FLAG >= 3)  Lz_ = sample_Lz; else Lz_ = 0; end;
		st = set_evolve_M_n_acq1Dgrad_t_CPP_struct(M,B_eff_rot,dta,dB0z,Ga_of_ta,sample_z_axis,RH_flag,...
		                                           acq_flag,Inhomo_Flag,Relax_Flag,T1,T2,context);
		[Mx,My,Mz,sig_real,sig_imag] = evolve_M_n_acq1Dgrad_t_CPP(st);
		M = [transpose(Mx),transpose(My),transpose(Mz)];
		Sig = transpose(sig_real) + i*transpose(sig_imag);
	else
		[M, Sig] = evolve_M_n_acq1Dgrad_t(M, B_eff_rot, dB0z, Ga_of_ta, dta, sample_z_axis, 1, Inhomo_Flag, 0);
	end;
else
	if (C_CODE)
		zGrad = ones(1,length(B_eff_rot(:,1)))*Ga(3);
		st = set_evolve_M_n_acq1Dgrad_t_CPP_struct(M,B_eff_rot,dta,dB0z,zGrad,sample_z_axis,RH_flag,...
		                                           acq_flag,Inhomo_Flag,Relax_Flag,T1,T2,context);
		[Mx,My,Mz,sig_real,sig_imag] = evolve_M_n_acq1Dgrad_t_CPP(st);
		M = [transpose(Mx),transpose(My),transpose(Mz)];
		Sig = transpose(sig_real) + i*transpose(sig_imag);
    else
    	[M, Sig] = evolve_M_n_acq(M, B_eff_rot, dB0z, Ga, dta, sample_z_axis, acq_flag, Inhomo_Flag, 0);
    end;
end;
plot_M_n_Phase(0,M_init,M,'M(z) post ACQ','Phase post ACQ',{'\phi(z)'},context);

% -------------------------------------------------------------------------
% Reconstruct signal
% -------------------------------------------------------------------------
M_final_z = abs(Sig / dZa);

if (~DesignPulse_Flag)
	R_rad   = 2*pi*R_of_te;                                                    % [rad/sec^2]
	Oi_rad  = 2*pi*Oi;                                                         % [rad/sec]
	alpha0 = (((omega_cs - Oi_rad)^2) ./ (2*R_rad)) - Tp*omega_cs;             % [rad]
	alpha1 = gammaHz*2*pi * Ge(3) * (((omega_cs - Oi_rad) ./ R_rad) - Tp);     % [rad/cm]
	alpha2 = ((gammaHz*2*pi)^2) * (Ge(3)^2) ./ (2*R_rad);                      % [rad/cm^2]
	
	za_axis = (gammaHz*2*pi * Ga(3) * ta - alpha1) / (2*alpha2);
end;

if (~DesignPulse_Flag && (DEBUG_FLAG >= 2))
	% Verify that the scanned range in z-axis is equal to the excited range
	z_chirp = O_of_te / (Ge(3)*gammaHz);

	figure; hold;
	plot((1:length(sample_z_axis)) /length(sample_z_axis) , sample_z_axis , 'b.-', ...
	     (1:length(z_chirp))       /length(z_chirp)       , z_chirp       , 'm.-', ...
	     (1:length(ta))            /length(ta)            , za_axis       , 'r.-');
	title('Comparison of Z-acquisition and Z-excitation');  xlabel('[none]');  ylabel('z [cm]');
	legend({'Sample z-axis','z-chirp','z-Acquisition'});
	grid; set_gca;
end;

figure; hold;
if (DesignPulse_Flag)
	plot(sample_z_axis, M_init(:,3),'k.-', hi_res_za, M_final_z(end:-1:1) ,'b.-');
else
	plot(sample_z_axis, M_init(:,3),'k.-', za_axis  , M_final_z(end:-1:1) ,'b.-');
end;
title('Initial & Reconstructed Magnetization');  xlabel('z [cm]');  ylabel('Magnetization');
axis([min(sample_z_axis),max(sample_z_axis),0,max(M_final_z)*1.2]); legend({'Initial M','Reconstructed M'});
set_gca; % set(gca,'FontSize',20);

% -------------------------------------------------------------------------
% Increase image resolution
% -------------------------------------------------------------------------
% Add noise
if (0)
	noise = transpose(rand(1,length(Sig)) + i*rand(1,length(Sig)));
	noise = noise - 0.5;
	noise = 0.05 * noise * max(abs(Sig));
	figure; plot(abs(Sig),'.-'); hold on;
	plot(abs(noise),'r.-');
	plot(abs(Sig+noise),'k.-');
	figure; plot(real(Sig),'.-'); hold on; plot(real(noise),'r.-');
	Sig2 = Sig + noise;
else
	Sig2 = Sig;
end;	

SR_M = increase_1D_image_resolution(Sig2(end:-1:1),alpha0,alpha1,alpha2,Ga(3),Ta,Lz,dZa,dZa/10);
SR_M = abs(SR_M) / max(abs(SR_M));
SR_z_axis = linspace(z_axis(1),z_axis(end),length(SR_M));

figure;
plot(sample_z_axis, M_init(:,3)         / max(M_init(:,3))        , 'k.-'); hold on;
plot(za_axis      , M_final_z(end:-1:1) / max(M_final_z(end:-1:1)), 'b.-');
plot(SR_z_axis    , SR_M                / max(SR_M)               , 'r.-');
legend({'S-Orig','S-Reconstructed','S-SR'});
title('Initial, Reconstructed & Super-Res Magnetization');  xlabel('z [cm]');  ylabel('Magnetization');

% -------------------------------------------------------
%  Create output file for Multi-SP post-processing stage
% -------------------------------------------------------
if (Nsp > 1)
	fn = sprintf('D:\\PhD\\Matlab\\Simulations\\Ver_4\\UF1D90\\MultSP_sim\\%3.0f.fid',FIDnum);
	fd = fopen(fn,'w');
	for idx = 1:length(Sig)
		fprintf(fd,'%5.5f %5.5f\n',real(Sig(idx)),imag(Sig(idx)));
	end;
	fclose(fd);
	
	fd = fopen(sprintf('D:\\PhD\\Matlab\\Simulations\\Ver_4\\UF1D90\\MultSP_sim\\set_%3.0f.m',FIDnum),'w');
	fprintf(fd,'sw   = %11.3f;       %% [kHz]   \n', 1e-3/dte     );
	fprintf(fd,'np   = %11.3f;       %%         \n', length(Sig)  );
	fprintf(fd,'Tp   = %11.3f;       %% [sec]   \n', Tp           );
	fprintf(fd,'Ta   = %11.3f;       %% [sec]   \n', Ta           );
	fprintf(fd,'Ge   = %11.3f;       %% [G/cm]  \n', Ge(3)        );
	fprintf(fd,'Ga   = %11.3f;       %% [G/cm]  \n', Ga(3)        );
	fprintf(fd,'Lz   = %11.3f;       %% [cm]    \n', Lz           );
	fprintf(fd,'rf_R = %11.3f;       %% [Hz/sec]\n', R_of_te      );
	fprintf(fd,'Nsp  = %11.3f;       %%         \n', Nsp          );
	fprintf(fd,'Zie  = %11.3f;       %% [cm]    \n', Zie          );
	fprintf(fd,'Zfe  = %11.3f;       %% [cm]    \n', Zfe          );
	fprintf(fd,'dZa  = %11.3f;       %% [cm]    \n', dZa          );
	
	fprintf(fd,'\n');
	fclose(fd);
end;

declare_end(sim_name); toc
return;

% ------------------------
% GE test code
% ------------------------
% if (0) % Gradient-Echo TEST - in order to view the phase, post excitation
%     declare_stage('Re-focus PI/2-Chirp pulse');
%     x_chirp_rot = zeros(1,length(chirp_rot));
%     y_chirp_rot = zeros(1,length(chirp_rot));
%     z_chirp_rot = zeros(1,length(chirp_rot));
%     B_eff_rot = [transpose(x_chirp_rot),transpose(y_chirp_rot),transpose(z_chirp_rot)] + (transpose(ones(1,length(chirp_rot))))*B_rot;
%     zGrad = -ones(1,length(B_eff_rot(:,1)))*Ge(3)/2;
%     st = set_evolve_M_n_acq1Dgrad_t_CPP_struct(M,B_eff_rot,dte,dB0z,zGrad,sample_z_axis,RH_flag,...
%                                                0,Inhomo_Flag,Relax_Flag,T1,T2,context);
%     [Mx,My,Mz,dummy1,dummy2] = evolve_M_n_acq1Dgrad_t_CPP(st);
%     M = [transpose(Mx),transpose(My),transpose(Mz)];
%     plot_M_n_Phase(0,M_init,M,'M(z) post refocusing pulse','Phase post refocusing pulse',{'\phi(z)'},context);
% 
%     declare_stage('Purge');
%     B_eff_rot = B_rot;
% 	st = set_evolve_M_CPP_struct(M,B_eff_rot,1.5E-3,dB0z,-20,sample_z_axis,RH_flag,...
% 	                             Inhomo_Flag,Relax_Flag,T1,T2,context);
% 	[Mx,My,Mz] = evolve_M_CPP(st);
% 	M = [transpose(Mx),transpose(My),transpose(Mz)];
%     plot_M_n_Phase(0,M_init,M,'M(z) post purge','Phase post purge',{'\phi(z)'},context);
%     
%     declare_stage('Acquiring signal');
%     Ta  = 3e-3;
%     dta = 2e-6;
%     ta  = 0:dta:(Ta-dta);
%     B_eff_rot = transpose(ones(1,length(ta))) * B_rot;
%     zGrad = 20 * ones(1,length(B_eff_rot(:,1)));
%     if (DEBUG_FLAG >= 3)  Lz_ = sample_Lz; else Lz_ = 0; end;
%     st = set_evolve_M_n_acq1Dgrad_t_CPP_struct(M,B_eff_rot,dta,dB0z,zGrad,sample_z_axis,RH_flag,...
%                                                1,Inhomo_Flag,Relax_Flag,T1,T2,context);
%     [Mx,My,Mz,sig_real,sig_imag] = evolve_M_n_acq1Dgrad_t_CPP(st);
%     M = [transpose(Mx),transpose(My),transpose(Mz)];
%     Sig = transpose(sig_real) + i*transpose(sig_imag);
%     plot_M_n_Phase(0,M_init,M,'M(z) post ACQ','Phase post ACQ',{'\phi(z)'},context);
%     
%     if (Inhomo_Fix_Flag)
%         save fix Sig;
%     else
%         save nofix Sig;
%     end;
% 
%     return;
% end;
