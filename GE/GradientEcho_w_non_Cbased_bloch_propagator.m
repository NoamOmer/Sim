% =================================
% Gradient-Echo
% =================================
clcl;
clear all;
simulation_name = 'Gradient-Echo Simulation';
declare_start(simulation_name);
compile_C_code();
context = 'ge_set_globals';
set_context;

% ------------------------------------------------
% Retrieve B0 inhomogeneity
% ------------------------------------------------
map_fn = [];
P_DNu0 = 0;
op     = 0;
dB0z   = zeros(1,length(sample_z_axis));

% -------------------------------------------------------------------------
% Create a shpaed sample
% -------------------------------------------------------------------------
[sample_z_axis,M_init] = UD1D_create_sample(context);
% return;

% -------------------------------------------------------------------------
% Excitation (non slice-selective hard pulse)
% -------------------------------------------------------------------------
declare_stage('Excitation');

pulse_op = 3;   % [1: default  |  2: SAR_test_f  |  3: load external pulse];

switch pulse_op
case 1
	pulse_prop.type    = Rect_P;
% 	pulse_prop.type    = Sinc_P;  pulse_prop.n_lobes = 4;  pulse_prop.Tp = 2*2.5e-3;    Gss = 1;  pulse_prop.calib   = 0;
	[B_RF_rot, B_RF_t] = gen_RF_pulse(pulse_prop, flip_a*pi/180, 0, dte, omega_CS, context);
	b_tmp = (B_RF_rot(:,1)+1i*B_RF_rot(:,2));
	figure; subplot(211); plot(B_RF_t,  abs(b_tmp),'.-'); title('amp');
			subplot(212); plot(B_RF_t,phase(b_tmp),'.-'); title('phs');
	Gss = 0;

case 2 % SAR test
	F = 1;
	pulse_prop.type    = Sinc_P;
	pulse_prop.n_lobes = 4;
	pulse_prop.Tp      = 2.5e-3/F;   % [sec]
	Gss                = 1*F;        % [G/cm]
	pulse_prop.calib   = 0;

	[B_RF_rot, B_RF_t] = gen_RF_pulse(pulse_prop, flip_a*pi/180, 0, dte, omega_CS, context);       % figure; plot(abs(B_RF_rot));
	b_tmp = (B_RF_rot(:,1)+1i*B_RF_rot(:,2))*gamma_T/(2*pi);
	pulse_flip_angle = sum(b_tmp)*dte*360;
	pulse_SAR        = sum((b_tmp')*(b_tmp))*dte;
	fprintf('Flip-angle = %3.3f [deg]\npulse_SAR = %3.3f [A.U.]\n',pulse_flip_angle,pulse_SAR);
	fprintf('%d deg excitation duration = %5.2f [us]\n',flip_a,(B_RF_t(end))*1E+6);

case 3 % external pulse
	fn = 'E:\01 Post\06 Projects\Siemens\04 extRF Pulse lib\SE2560A180.SE180_12A2_2.pta';
	Tp = 2560e-3;
	[b1, N] = read_siemens_PTA_pulse(fn);
	dte = Tp/N;
	B_RF_t = 0:dte:Tp-dte;
	b1 = rf_calib_by_pulse_integral(b1,dte,flip_a*2);
	
	B_RF_rot(:,1) = real(b1);
	B_RF_rot(:,2) = imag(b1);
	B_RF_rot(:,3) = zeros(length(b1),1);

	Gss = 0;
	
	figure;
	subplot(211); plot(  abs(b1),'.-'); title('Pulse Amplitude');
	subplot(212); plot(angle(b1),'.-'); title('Pulse phase [rad]');
end;


B_eff_rot = B_RF_rot + (transpose(ones(1,length(B_RF_t)))) * B_rot;
if (C_CODE)
 	Exc_Gradient = Gss*ones(1,length(B_eff_rot(:,1)));
  	st = set_evolve_M_n_acq1Dgrad_t_CPP_struct(M_init,B_eff_rot,dte,dB0z,Exc_Gradient,sample_z_axis,...
  	                                           RH_flag,0,Inhomo_Flag,Relax_Flag,T1,T2,context);
 	[Mx,My,Mz,dummy1,dummy2] = evolve_M_n_acq1Dgrad_t_CPP(st);
 	M = [transpose(Mx),transpose(My),transpose(Mz)];
else
    Exc_Gradient = [0,0,0];
    M = evolve_M_n_acq(M_init,B_eff_rot,dB0z,Exc_Gradient,dt,sample_z_axis,0,Inhomo_Flag,0);
end;
plot_M_n_Phase(0,M_init,M,'M(z) after pi/2 excitation pulse','Phase after pi/2 excitation',{'\phi(z)'},context);


% -------------------------------------------------------------------------
% Pre-acquisition refocusing gradient
% -------------------------------------------------------------------------
declare_stage('Pre-acq refocusing gradient'); tic
B_eff_rot = B_rot;
if (C_CODE)
	st = set_evolve_M_CPP_struct(M,B_eff_rot,Tpr,dB0z,Gpr,sample_z_axis,RH_flag,...
	                             Inhomo_Flag,Relax_Flag,T1,T2,context);
	[Mx,My,Mz] = evolve_M_CPP(st);
 	M = [transpose(Mx),transpose(My),transpose(Mz)];
else
    M = evolve_M(M,B_eff_rot,dB0z,[0,0,Gpr],Tpr,sample_z_axis,Inhomo_Flag);
end;
% plot_M_n_Phase(0,M_init,M,'M(z) after Pre-acq refocusing gradient','Phase after Pre-acq refocusing gradient',{'\phi(z)'},context);


% -------------------------------------------------------------------------
% Acquisition
% -------------------------------------------------------------------------
declare_stage('Acquisition');
B_eff_rot = transpose(ones(1,length(ta))) * B_rot;
acq_flag = 1;
if (C_CODE)
    acq_grad = Ga * ones(1,length(ta));
    st = set_evolve_M_n_acq1Dgrad_t_CPP_struct(M,B_eff_rot,dta,dB0z,acq_grad,sample_z_axis,RH_flag,acq_flag,...
                                               Inhomo_Flag,Relax_Flag,T1,T2,context);
    [Mx,My,Mz,sig_real,sig_imag] = evolve_M_n_acq1Dgrad_t_CPP(st);
    M = [transpose(Mx),transpose(My),transpose(Mz)];
    Sig = transpose(sig_real) + 1i*transpose(sig_imag);
else
    Exc_Gradient = [0,0,1]*Ga;
    [M,Sig] = evolve_M_n_acq(M,B_eff_rot,dB0z,Exc_Gradient,dta,sample_z_axis,1,Inhomo_Flag,0);
    sig_real = real(Sig);
    sig_imag = imag(Sig);
end;

% -------------------------------------------------------------------------
% Add noise
% -------------------------------------------------------------------------
FT_E                = [];
orig_sig            = transpose(Sig);
reset_random_seed_f = 1;
interp_flag         = 1;
% % % % % for idx = 1:256, fprintf('%1.0f\n',idx);
% % % % % if (1)
% % % % % % 	declare_stage('Adding noise');
% % % % % 	noise_std           = 1;             % [%]
% % % % % 	signal_lvl          = 1;
% % % % % 	complex_noise_f     = 1;
% % % % % 	Sig = add_random_noise_to_1D_sig(orig_sig,fb,noise_std,signal_lvl,complex_noise_f,reset_random_seed_f,0);
% % % % % 	reset_random_seed_f = 0;
% % % % % end;

% -------------------------------------------------------------------------
% Post processing - ZF, FFT & display
% -------------------------------------------------------------------------
% declare_stage('Post-Processing');

if (interp_flag)
	Sig2 = zero_pad_1D(zero_pad_1D(zero_pad_1D(zero_pad_1D(Sig))));
else
	Sig2 = Sig;
end;
za_axis = linspace(za_axis(1),za_axis(end),length(Sig2));

% final_image = fftshift(fft(fftshift(Sig2)))*sqrt(N_acq_pts);

if (isempty(FT_E)), FT_E = calc_transformation_SNR_change__calcExplicitFTmatrix(length(Sig2)); end;
final_image = flipud(fftshift(ctranspose(FT_E) * transpose(Sig2)));

% final_image = abs(final_image); final_image = final_image / max(final_image);

% % % % % im_series{idx}  = final_image;
% % % % % sig_series{idx} = Sig;
% % % % % end;
% % % % % 
% % % % % save('8_2_1_im_series_ge','sig_series','im_series');

final_image = abs(final_image);
% final_image = final_image / max(final_image);
M0 = M_init(:,3) * max(final_image) / max(M_init(:,3));

figure;
if (1)
	plot(za_axis      ,final_image,'b.-');
% 	axis([-6 6 2.81 2.83]); zoom;
	title(sprintf('Final image (N_{pts} = %3.0f;  Ga = %3.4f [G/cm];  Ta = %3.1f [ms])',N_acq_pts,Ga,Ta*1e3));	xlabel('[cm]');
else
	subplot(211); hold on;
				  plot(abs(Sig),'k.-');
% 				  plot(real(Sig),'g.-');
% 				  plot(imag(Sig),'b.-');
				  title(sprintf('Acquired signal')); legend({'abs','real','imag'});
	subplot(212); hold on;
				  plot(sample_z_axis,M0         ,'k-' );
				  plot(za_axis      ,final_image,'b.-');
				  title(sprintf('Final image (N_{pts} = %3.0f;  Ga = %3.4f [G/cm])',N_acq_pts,Ga));  xlabel('[cm]');
end;

declare_end(simulation_name);

return;

% plot the PSF
z   = sample_z_axis;
Lz  = 2*Lhalf;
dza = Lz/N_acq_pts;
PSF = (exp(1i*pi*z/Lz) .* sin(pi*z/dza)) ./ (sin(pi*z/Lz));
acq_grid = linspace(za_axis(1),za_axis(end),N_acq_pts);

figure; hold on;
plot(z,abs(PSF),'b-');
plot(acq_grid,zeros(1,N_acq_pts),'m^-');
legend({'PSF','acquisition grid'});
