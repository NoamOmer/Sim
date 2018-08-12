% 1D Spatiotemporally-Encoded MRI Simulation
function UF1D()
clear all; pack;
% clcl;
tic
sim_name = '1D SPEN MRI';
declare_start(sim_name);  tic
compile_C_code();
context = 'UF_1D_set_globals';
set_context; %set_globals;

% % -------------------------------------------------------------------------
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
[~,P_DNu0,~] = UD1D_retrieve_inhomogeneity(context);

% -------------------------------------------------------------------------
% Create a shaped sample
% -------------------------------------------------------------------------
[sample_z_axis,M_init] = UD1D_create_sample(context);
% return;

% Now that the z_axis is arranged -- define the inhomogeneity field
% B_lab = B0_lab + [0, 0, (omega_CS /gamma_T)];
% B_rot = B_lab  - [0, 0, (omega_rot/gamma_T)];
B_rot = [0,0,0];
dOmega_0 = polyval(P_DNu0,sample_z_axis - 0*sample_z_axis(length(sample_z_axis)/8));      % [Hz]
dB0z = (dOmega_0/gammaHz)*1E-4;                                                           % [T]
if (sum(P_DNu0) > 0)
figure; plot(sample_z_axis*10,dB0z*gamma_Hz_T*1e-3,'.-'); title('dB_0'); ylabel('kHz'); xlabel('[mm]'); grid;
end;

% -------------------------------------------------------------------------
% Design & Apply PI/2-Chirp pulse
% -------------------------------------------------------------------------
[chirp_rot,phi_chirp_rot,alpha0,alpha1,alpha2] = UD1D_design_pulse(context);
chirp_rot = chirp_rot * exp(1i*exc_phase);
SAR_test_f = 0;
if (SAR_test_f)
pulse_flip_angle = sum(chirp_rot)*dte*360;
pulse_SAR        = sum((chirp_rot')*(chirp_rot))*dte;
fprintf('Flip-angle = %3.3f [deg]\npulse_SAR = %3.3f [A.U.]\n',pulse_flip_angle,pulse_SAR);
end;

if (0)
	[chirp_rot,phi_chirp_rot] = load_external_RF('E:\01 Post\06 Projects\Olfaction\2011_12_29_3T_wig_rand_etc\chirp_RF_Ge8.txt');
end;

% verify_inhom_fix_1D(context);

x_chirp_rot = real(chirp_rot);
y_chirp_rot = imag(chirp_rot);
z_chirp_rot = zeros(1,length(chirp_rot));

if (DEBUG_FLAG >= 3)
	te_plot = linspace(te(1), te(end), length(x_chirp_rot));
	figure; hold on;
	subplot(311); plot(te_plot*1e+3, x_chirp_rot, 'k.-', te_plot*1e+3, y_chirp_rot, 'b.-');
	title('Chirp pulse in rotating frame');      xlabel('Time [ms]');  ylabel('Magnetic field'); 
	legend({'x-coordinate','y-coordinate'});     set_gca;
	subplot(312); plot(te_plot*1e+3, abs(chirp_rot), 'b.-');
	title('Chirp amplitude in rotating frame');  xlabel('Time [ms]');  ylabel('Amp [T]'); set_gca;
	subplot(313); plot(te_plot*1e+3, phi_chirp_rot, 'r.-');
	title('Chirp angle in rotating frame');      xlabel('Time [ms]');  ylabel('Angle [rad]'); set_gca;
end;

declare_stage('Apply PI/2-Chirp pulse');
B_eff_rot = [transpose(x_chirp_rot),transpose(y_chirp_rot),transpose(z_chirp_rot)] + (transpose(ones(1,length(chirp_rot))))*B_rot;
if (C_CODE)
	zGrad = ones(1,length(B_eff_rot(:,1)))*Ge(3);
	st = set_evolve_M_n_acq1Dgrad_t_CPP_struct(M_init,B_eff_rot,dte,dB0z,zGrad,sample_z_axis,RH_flag,...
	                                           0,Inhomo_Flag,Relax_Flag,T1,T2,context);
	[Mx,My,Mz,~,~] = evolve_M_n_acq1Dgrad_t_CPP(st);
	M = [transpose(Mx),transpose(My),transpose(Mz)];
else
	[M,~] = evolve_M_n_acq(M_init, B_eff_rot, dB0z, Ge, dte, sample_z_axis, 0, Inhomo_Flag, 0);
end;
plot_M_n_Phase(0  ,sample_z_axis,M_init,M,'M(z) post PI/2 Chirp pulse','Phase post PI/2 Chirp pulse',{'\phi(z)'},context);

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

	zGrad = ones(1,length(B_eff_rot(:,1))) * (-Ge(3)*Ge_AUGMENT_F);
	st = set_evolve_M_n_acq1Dgrad_t_CPP_struct(M,B_eff_rot,dt_refocus,dB0z,zGrad,sample_z_axis,RH_flag,...
											   acq_flag,Inhomo_Flag,Relax_Flag,T1,T2,context);
	[Mx,My,Mz,~,~] = evolve_M_n_acq1Dgrad_t_CPP(st);
	M = [transpose(Mx),transpose(My),transpose(Mz)];
	plot_M_n_Phase(0,sample_z_axis,M_init,M,'M(z) post SLR refocusing','Phase post SLR refocusing',{'\phi(z)'},context);
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
	[Mx,My,Mz,~,~] = evolve_M_n_acq1Dgrad_t_CPP(st);
	M = [transpose(Mx),transpose(My),transpose(Mz)];
	plot_M_n_Phase(0,sample_z_axis,M_init,M,'M(z) Purge 1st part','Phase Purge 1st part',{'\phi(z)'},context);

	% Reverse the acquisition gradient
	Ga_of_ta = -Ga_of_ta;
end;

% -------------------------------------------------------------------------
% Spin-Echo pulse
% -------------------------------------------------------------------------
if (SE_flag)
% 	dt_crush   = dte;
% 	nCrush_pts = Tcrush/dt_crush;
% 	x_chirp_rot = zeros(1,nCrush_pts);
% 	y_chirp_rot = zeros(1,nCrush_pts);
% 	z_chirp_rot = zeros(1,nCrush_pts);
% 	B_crush_rot = [transpose(x_chirp_rot),transpose(y_chirp_rot),transpose(z_chirp_rot)];
% 	Gcrush = Gcrush * ones(1,length(B_crush_rot(:,1)));

	% Pre PI crusher
	if (Tcrush ~= 0)
		declare_stage('Pre PI Pulse Crusher');
		B_crush_rot = B_rot;
		st = set_evolve_M_CPP_struct(M,B_crush_rot,Tcrush,dB0z,Gcrush,sample_z_axis,RH_flag,Inhomo_Flag,Relax_Flag,T1,T2,context);
		[Mx,My,Mz] = evolve_M_CPP(st);		
% 		st = set_evolve_M_n_acq1Dgrad_t_CPP_struct(M,B_crush_rot,dt_crush,dB0z,Gcrush,sample_z_axis,RH_flag,...
% 												   0,Inhomo_Flag,Relax_Flag,T1,T2,context);
% 		[Mx,My,Mz,~,~] = evolve_M_n_acq1Dgrad_t_CPP(st);
		M = [transpose(Mx),transpose(My),transpose(Mz)];
	end;
	
	% PI pulse
	refocus_pulse_op = 1;
	switch (refocus_pulse_op)
	case 1
		declare_stage('PI Pulse (non slice--selective hard pulse)');
		dte_refoc        = 10e-6;
		pulse_prop.type  = Rect_P;
		pulse_prop.calib = 1;
		[B_RF_rot, B_RF_t] = gen_RF_pulse(pulse_prop, pi, refoc_phase, dte_refoc, omega_CS, context);
		fprintf('PI Pulse duration = %d [us]\n',(B_RF_t(end))*1E+6);
		B_eff_rot = B_RF_rot + (transpose(ones(1,length(B_RF_t)))) * B_rot;
	case 2
		dte_refoc = Trefoc/N_refoc;
		pulse_prop.type = Sinc_P;    pulse_prop.n_lobes = n_lobes_refoc;    pulse_prop.Tp = Trefoc;    pulse_prop.calib = 1;
		[B_RF_rot, B_RF_t] = gen_RF_pulse(pulse_prop, (refoc_angle)*pi/180, refoc_phase, dte_refoc, omega_CS, context);
		B_eff_rot = B_RF_rot + (transpose(ones(1,length(B_RF_t)))) * B_rot;	
	end;
	acq_flag = 0;

	zGrad = zeros(1,length(B_RF_rot(:,1)));
	st = set_evolve_M_n_acq1Dgrad_t_CPP_struct(M,B_eff_rot,dte_refoc,dB0z,zGrad,sample_z_axis,RH_flag,...
											   acq_flag,Inhomo_Flag,Relax_Flag,T1,T2,context);
	[Mx,My,Mz,~,~] = evolve_M_n_acq1Dgrad_t_CPP(st);
	M = [transpose(Mx),transpose(My),transpose(Mz)];

	% Post PI crusher
	if (Tcrush ~= 0)
		declare_stage('Post PI Pulse Crusher');
		st = set_evolve_M_CPP_struct(M,B_crush_rot,Tcrush,dB0z,Gcrush,sample_z_axis,RH_flag,Inhomo_Flag,Relax_Flag,T1,T2,context);
		[Mx,My,Mz] = evolve_M_CPP(st);
% 		st = set_evolve_M_n_acq1Dgrad_t_CPP_struct(M,B_crush_rot,dt_crush,dB0z,Gcrush,sample_z_axis,RH_flag,...
% 												   0,Inhomo_Flag,Relax_Flag,T1,T2,context);
% 		[Mx,My,Mz,~,~] = evolve_M_n_acq1Dgrad_t_CPP(st);
		M = [transpose(Mx),transpose(My),transpose(Mz)];
	end;
	
	plot_M_n_Phase(0,sample_z_axis,M_init,M,'M(z) after PI pulse','Phase after PI pulse',{'\phi(z)'},context);
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
	[Mx,My,Mz,~,~] = evolve_M_n_acq1Dgrad_t_CPP(st);
	M = [transpose(Mx),transpose(My),transpose(Mz)];
	plot_M_n_Phase(0,sample_z_axis,M_init,M,'M(z) Purge 2nd part','Phase Purge 2nd part',{'\phi(z)'},context);
end;

% -------------------------------------------------------------------------
% Acquire signal
% -------------------------------------------------------------------------
declare_stage('Acquiring signal');
B_eff_rot = transpose(ones(1,length(ta))) * B_rot;
acq_flag = 1;
if (DesignPulse_Flag)
	if (C_CODE)
		st = set_evolve_M_n_acq1Dgrad_t_CPP_struct(M,B_eff_rot,dta,dB0z,Ga_of_ta,sample_z_axis,RH_flag,...
		                                           acq_flag,Inhomo_Flag,Relax_Flag,T1,T2,context);
		[Mx,My,Mz,sig_real,sig_imag] = evolve_M_n_acq1Dgrad_t_CPP(st);
		M = [transpose(Mx),transpose(My),transpose(Mz)];
		Sig = transpose(sig_real) + 1i*transpose(sig_imag);
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
		Sig = transpose(sig_real) + 1i*transpose(sig_imag);
    else
    	[M, Sig] = evolve_M_n_acq(M, B_eff_rot, dB0z, Ga, dta, sample_z_axis, acq_flag, Inhomo_Flag, 0);
    end;
end;
Sig = transpose(Sig);
plot_M_n_Phase(0,sample_z_axis,M_init,M,'M(z) post ACQ','Phase post ACQ',{'\phi(z)'},context);
% -------------------------------------------------------------------------
% Reconstruct signal
% -------------------------------------------------------------------------
if (~DesignPulse_Flag)
	R_rad   = 2*pi*rate;                                                       % [rad/sec^2]
	Oi_rad  = 2*pi*Oi_Hz;                                                      % [rad/sec]
	alpha0 = (((omega_CS - Oi_rad)^2) ./ (2*R_rad)) - Tp*omega_CS;             % [rad]
	alpha1 = gammaHz*2*pi * Ge(3) * (((omega_CS - Oi_rad) ./ R_rad) - Tp);     % [rad/cm]
	alpha2 = ((gammaHz*2*pi)^2) * (Ge(3)^2) ./ (2*R_rad);                      % [rad/cm^2]
	
	za_axis = (gammaHz*2*pi * Ga(3) * ta - alpha1) / (2*alpha2);
end;

if (DesignPulse_Flag)
	plot(sample_z_axis*10, M_init(:,3),'k-', hi_res_za*10, abs(Sig(end:-1:1)) ,'b-');
end;

% -------------------------------------------------------------------------
% Increase image resolution
% -------------------------------------------------------------------------
% declare_stage('Increasing signal resolution');

% SR_Factor   = F_SR; this is not important - we just ask the SR recon to produce an image with N_acq_pts
% SR_Factor   = SR_Factor*256/(length(Sig));

% GePE        = Ge(3);    NSE         = length(Sig);
% TaPE        = dta;      GaPE        = -Ga(3);
SRExp       = 0;%14;
SRESP       = 2;
SR_HRF      = 15;
SR_resp_win = 1;                          % [1: PSF               |  2: Gauss]
SR_ones_f   = 0;
% SR_cs.exp   = 0;
% SR_cs.sft   = 0;        SR_cs.c2    = 0;
% SR_cs.ESP   = 2;        SR_cs.win   = 1;  % [1: SPEN based gauss  |  2: simple smooth filter  |  3: orig exponential filter]

ax.fh_PSF      = 0;%figure;
ax.ax_A        = 0;
ax.ax_W_f      = 0;
ax.ax_spect    = 0;
ax.ax_A_Atag   = 0;
ax.ax_AW_AWtag = 0;
ax.ax_EIG      = 0;

% == Add noise ==
SPEN_B              = [];
SPEN_E              = [];
reset_random_seed_f = 1;
interp_flag         = 1;

% % % % % for idx = 1:256, fprintf('%1.0f\n',idx);
% % % % % if (1)
% % % % % % 	declare_stage('Adding noise');
% % % % % 	noise_std           = 1;             % [%]
% % % % % 	signal_lvl          = 1;
% % % % % 	complex_noise_f     = 1;
% % % % % 	[Sig,~] = add_random_noise_to_1D_sig(orig_sig,fb,noise_std,signal_lvl,complex_noise_f,reset_random_seed_f,0);
% % % % % 	reset_random_seed_f = 0;
% % % % % end;

% disp(sprintf('Trying to increase image resolution %3.3f [mm] by a factor of %3.1f',dZa,SR_Factor));

% [nbe] replace manual SR with standard one
% SR_Im = increase_1D_image_resolution(Sig(end:-1:1),alpha0,alpha1,alpha2,Ga(3),Ta,Lz_imag,dZa,dZa/SR_Factor,SRExp,SRESP,SR_HRF);

% [SR_Im,~]=increase_2D_Hybrid_image_resolution5(Sig,...
%     Zie,Zfe,Lz_imag,dZa,dZa/SR_Factor,alpha0/(2*pi),alpha1/(2*pi),alpha2/(2*pi),GePE,Tp,NSE,...
%     TaPE,Ta,GaPE,SE_flag,SRExp,SRESP,0,1,SR_HRF,1,0,SR_cs,ax,0,0,SR_resp_win,0,[]);
% SR_Im = SR_Im*N_acq_pts;

% SRExp       = 0;
% SR_resp_win = 1;                          % [1: PSF               |  2: Gauss]
% SR_ones_f   = 0;
if (interp_flag),  req_N_acq_pts = 512;         nonSR_Im = interp1(1:length(Sig),abs(Sig),linspace(1,length(Sig),req_N_acq_pts));
else               req_N_acq_pts = length(Sig); nonSR_Im = Sig;  end;

% if (isempty(SPEN_B) || isempty(SPEN_E))
	plot_f = 0;
	[SPEN_B SPEN_E] = calc_transformation_SNR_change__calcExplicitSPENmatrix(req_N_acq_pts,length(Sig),Lz_imag,Tp,Ta,Ge(3),Ga(3),...
	                                                                         SE_flag,SRExp,SRESP,SR_HRF,SR_ones_f,gammaHz,...
																			 SR_resp_win,plot_f);
% end;
SR_Im = transpose(SPEN_B*transpose(Sig));
% figure; 
% subplot(211); plot(abs(SR_Im),'.-');                                          title('SR');
% subplot(212); plot(sample_z_axis*10, M_init0,'k-', za_axis*10, Im_FFT,'b-');  title('fft'); legend({'M_0','Image'});

% % % % % sig_series{idx} = reverse(Sig);
% % % % % im_series {idx} = SR_Im;
% % % % % end; % noise forloop
% % % % % save(sprintf('8_2_1_im_series_spen_Fsr_%2.2f.mat',F_SR),'sig_series','im_series');

nonSR_Im = abs(nonSR_Im(end:-1:1));
SR_SigPh = phase(SR_Im);
SR_Im    = abs(SR_Im);

% ----------------
%  FT processing
% ----------------
if (interp_flag)
% 	Sig2 = zero_pad_1D(zero_pad_1D(zero_pad_1D(zero_pad_1D(zero_pad_1D(Sig)))));
	Sig2 = zero_pad_1D(zero_pad_1D(zero_pad_1D(Sig)));
else
	Sig2 = Sig;
end;

FT_proc = 1;
if (FT_proc)
% 	R_HzSec    = gammaHz*Ge(3)*Lz/Tp;
% 	alpha2     = - ((gammaHz*Ge(3))^2) / (2*R_HzSec);           % [1/cm^2]
% 	Sig2 = Sig2 .* transpose(exp(i*alpha2*(za_axis.^2)));
	Im_FFT = abs(fftshift(fft(fftshift(Sig2))));
	za_axis = linspace(Zfa,Zia,length(Im_FFT));
% 	figure; hold;
% 	plot(sample_z_axis*10, M_init0,'k-', za_axis*10, Im_FFT,'b-');
end;

% ----------------
%  Normalization
% ----------------
normalize_Im_flag = 0;
if (normalize_Im_flag)
Im_FFT   = abs(Im_FFT) / max(abs(Im_FFT));
nonSR_Im = nonSR_Im    / max(nonSR_Im);
SR_Im    = SR_Im       / max(SR_Im);
else
% M0_nonSR = M_init(:,3);
% M0_SR    = M_init(:,3);
end;
M0_nonSR = M_init(:,3) * max(nonSR_Im) / (max(M_init(:,3)));
M0_SR    = M_init(:,3) * max(SR_Im)    / (max(M_init(:,3)));
% M0_fft   = M_init(:,3) * max(Im_FFT)   / (max(M_init(:,3)));

% ----------------
%  plot
% ----------------
% figure;
% za_axis2 = linspace(Zfa,Zia,length(nonSR_Im));
% subplot(211); plot(za_axis2*10,abs(nonSR_Im),' -b');  hold on;
% 			  plot(sample_z_axis*10,M0_nonSR,' -k');  title(['nonSR_Im (\Deltax = ' sprintf('%3.2f [mm])',dZa*10)]);  grid;
% 			  xlabel('x-axis [mm]'); legend({'non SR','M_0'});
% za_axis2 = linspace(Zfa,Zia,length(SR_Im));
% % 	[~,center_idx] = max(SR_Im);
% % 	center_integral = sum(SR_Im(center_idx-1:center_idx+1)) / 3;
% % 	center_integral = SR_Im(center_idx);
% subplot(212); plot(za_axis2*10,SR_Im       ,' -b');   hold on;
% 			  plot(sample_z_axis*10,M0_SR  ,' -k');
%               plot(za_axis*10,Im_FFT       ,' -r');
% 			  title('SR'); xlabel('x-axis [mm]'); grid; legend({'SR','M_0','FFT'});
% % 			  axis([-Lhalf*10/2 +Lhalf*10/2 0 1.1]);
% % subplot(223); plot(sample_z_axis*10, M0_fft,'-k', za_axis*10, abs(Im_FFT),'-b');  title('fft'); legend({'M_0','FFT'});

% ================ for PAPER: compare Nyquist vs. non-Nyquist
figure;
za_axis2 = linspace(Zfa,Zia,length(SR_Im));
plot(sample_z_axis ,M0_SR ,' -k'); hold on;
plot(za_axis2      ,SR_Im ,'.-b'); % set(gca,'xticklabel',[],'yticklabel',[]);
axis([-20    20     0    35]);

% % % % % ================ for PAPER: proove SPA
% % % % figure;
% % % % za_axis2 = linspace(Zfa,Zia,length(SR_Im));
% % % % plot(sample_z_axis*10 ,M0_nonSR ,' -k'); hold on;
% % % % plot(za_axis2*10      ,nonSR_Im ,'.-b'); set(gca,'xticklabel',[],'yticklabel',[]);
% % % % title(sprintf('SPEN - non SR (req res = %3.3f [mm])',required_res*10));


% ================ for PAPER: compare Nyquist vs. non-Nyquist
% % % % figure;
% % % % za_axis2 = linspace(Zfa,Zia,length(SR_Im));
% % % % plot(sample_z_axis*10 ,M0_SR ,' -k'); hold on;
% % % % plot(za_axis2*10      ,SR_Im ,'.-b'); set(gca,'xticklabel',[],'yticklabel',[]);
% % % % title('@ Nyquist (dta x1.0) - SPEN - SR wOnes');
% % % % % title('non Nyquist (dta x2.5) - SPEN - SR');
% % % % 
% % % % figure;
% % % % plot(sample_z_axis*10 ,M0_SR ,' -k'); hold on;
% % % % plot(za_axis*10       ,Im_FFT,'.-b'); set(gca,'xticklabel',[],'yticklabel',[]);
% % % % title('@ Nyquist (dta x1.0) - SPEN - FFT');
% % % % % title('non Nyquist (dta x1.0) - SPEN - FFT');


% -------------------------------------------------------
%  Create output file for Multi-SP post-processing stage
% -------------------------------------------------------
if (Nsp > 1)
	fn = sprintf('D:\\PhD\\Matlab\\Simulations\\Ver_4\\UF1D90\\MultSP_sim\\%3.0f.fid',FIDnum);
	fd = fopen(fn,'w');
	parfor idx = 1:length(Sig)
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
	fprintf(fd,'Lz   = %11.3f;       %% [cm]    \n', Lz_imag      );
	fprintf(fd,'rf_R = %11.3f;       %% [Hz/sec]\n', R_of_te      );
	fprintf(fd,'Nsp  = %11.3f;       %%         \n', Nsp          );
	fprintf(fd,'Zie  = %11.3f;       %% [cm]    \n', Zie          );
	fprintf(fd,'Zfe  = %11.3f;       %% [cm]    \n', Zfe          );
	fprintf(fd,'dZa  = %11.3f;       %% [cm]    \n', dZa          );
	
	fprintf(fd,'\n');
	fclose(fd);
end;

declare_end(sim_name); toc

toc
return;


% ================================================================================================
% ================================================================================================


function [chirp_rot,phi_chirp_rot] = load_external_RF(fn)
len = 512;
a=load(fn);
loaded_chirp_rot     = transpose(max(abs(chirp_rot)) * (a(1:len,1) / max(abs(a(1:len,1)))));
loaded_phi_chirp_rot = transpose(a(1:len,2));

% figure;plot(loaded_chirp_rot,'k.-');
% figure;plot(unwrap(loaded_phi_chirp_rot),'k.-');

figure;
subplot(211); plot(abs(chirp_rot),'k.-'); hold on;
subplot(212); plot(phi_chirp_rot ,'k.-'); hold on;
chirp_rot = loaded_chirp_rot .* exp(1i*loaded_phi_chirp_rot);
phi_chirp_rot = unwrap(loaded_phi_chirp_rot);
subplot(211); plot(abs(chirp_rot),'r.-'); legend({'sim','poet'});
subplot(212); plot(phi_chirp_rot ,'r.-'); legend({'sim','poet'});

return;


% ================================================================================================
% ================================================================================================



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
%     plot_M_n_Phase(0,sample_z_axis,M_init,M,'M(z) post refocusing pulse','Phase post refocusing pulse',{'\phi(z)'},context);
% 
%     declare_stage('Purge');
%     B_eff_rot = B_rot;
% 	st = set_evolve_M_CPP_struct(M,B_eff_rot,1.5E-3,dB0z,-20,sample_z_axis,RH_flag,...
% 	                             Inhomo_Flag,Relax_Flag,T1,T2,context);
% 	[Mx,My,Mz] = evolve_M_CPP(st);
% 	M = [transpose(Mx),transpose(My),transpose(Mz)];
%     plot_M_n_Phase(0,sample_z_axis,M_init,M,'M(z) post purge','Phase post purge',{'\phi(z)'},context);
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
%     plot_M_n_Phase(0,sample_z_axis,M_init,M,'M(z) post ACQ','Phase post ACQ',{'\phi(z)'},context);
%     
%     if (Inhomo_Fix_Flag)
%         save fix Sig;
%     else
%         save nofix Sig;
%     end;
% 
%     return;
% end;
