% 1D Ultrafast MRI Simulation
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
	disp(sprintf('\n  Dk = [%3.1f ... %3.1f]\n  k0 = %3.1f\n',Dk/2,+Dk/2,k0_tag));
end;

% -------------------------------------------------------------------------
% Retrieve the B0 inhomogeneity and set the spatial axis
% -------------------------------------------------------------------------
declare_stage('Retrieve the B0 inhomogeneity and set the spatial axis');
fh = 0;
if (Map_Mult_Factor ~= 0)
	if (Seq1D)
		[map_fn,P_DNu0,z_axis,tof,op,mean_std] = Inhomo_Zmap;
	else
		map_fn = '31mar08_1_MATOrient_RP90_SF2_DSF1_PO5.mat';%uigetfile('*.mat','Select an Oriented Phase-Map File');
		load(map_fn);
		[z_axis,P_DNu0,Porder] = calc_oriented_1D_inhomo(map_mat, map_mask, x_axis, y_axis, rot_phi, 0, 0, gPolyOrder);
	end;
	if (DEBUG_FLAG >= 1)
		fh = figure; plot(z_axis,polyval(P_DNu0,z_axis)); hold on;
		title('Field inhomogeneity'); xlabel('z-axis [cm]'); ylabel('\Delta\Nu [Hz]'); set_gca;
	end;
	% Inhomogeneous field tests
	P_DNu0 = Map_Mult_Factor*P_DNu0;
	plot(z_axis,polyval(P_DNu0,z_axis),'c');
else
	P_DNu0 = 0;
end;

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
% Create a new z-axis for the the sample
sample_z_axis = (sample_Zi : dz : sample_Zf+dz);

% % Triple Gaussian
% loc = 400;
% sample = exp(-2E+1*((sample_z_axis+loc*dz - 50*dz).^2)) + ...
% 	     exp(-2E+1*((sample_z_axis-loc*dz - 50*dz).^2)) + ...
% 		 exp(-2E+1*((sample_z_axis- 50*dz - 50*dz).^2));      % Smooth shape (gaussian) decreases the effect of Gibbs phenomenon

% Gaussian
sample = exp(-70E+21*((sample_z_axis *dz).^8));

% % Square
% sample = [zeros(1,750), ones(1,length(sample_z_axis)-1500), zeros(1,750)];   % This shape is convenient for examining the excitation efficiency

% % Triangle
% F = 500;
% sample = [zeros(1,F), 1:(length(z_axis)/2-F), (length(z_axis)/2-F):-1:1, zeros(1,F)]; sample = sample/max(sample);
% if (length(sample) < length(z_axis))   sample(end+1) = sample(end);   end;

% % Ones
% sample = ones(1,length(sample_z_axis));

sample = sample ./ max(sample);
sample = sample + 1E-10;
M_init = transpose(sample) * M0;
if (DEBUG_FLAG >= 0)
	figure; plot(sample_z_axis,M_init(:,3),'.-'); title('Sample(z)'); xlabel('z-axis [cm]'); ylabel('M'); set_gca;
end;

% Now that the z_axis is arranged -- define the inhomogeneity
B_lab = B0_lab + [0, 0, (omega_cs /gamma_T)];
B_rot = B_lab  - [0, 0, (omega_rot/gamma_T)];
dOmega_0 = polyval(P_DNu0,sample_z_axis);      % [Hz]
dB0z = (dOmega_0/gammaHz)*1E-4;                % [T]

% -------------------------------------------------------------------------
% Design & Apply PI/2-Chirp pulse
% -------------------------------------------------------------------------
declare_stage('Design PI/2-Chirp pulse');
if (DesignPulse_Flag)
	if (Inhomo_Fix_Flag == 0)
		P_DNu0_ = P_DNu0;
		P_DNu0  = P_DNu0*0;
	end;
	if (RH_flag)
		if (SE_flag)
			[te,dte,O_of_te,Ge,ta,Ga_of_ta,R_of_te,chirp_rot,phi_chirp_rot]=design_inhom_fix_1DUF_RHR_Ver9_1(context);
		else
			[te,dte,O_of_te,Ge,ta,Ga_of_ta,R_of_te,chirp_rot,phi_chirp_rot]=design_inhom_fix_1DUF_RHR_Ver9(context);
		end;
	else
		[te,dte,O_of_te,Ge,ta,Ga_of_ta,R_of_te,chirp_rot,phi_chirp_rot] = design_inhom_fix_1DUF(context);
	end;
	if (Inhomo_Fix_Flag == 0)
		P_DNu0 = P_DNu0_;
	end;
    DOkHz = 1E-3*abs(O_of_te(end) - O_of_te(1));
    disp(sprintf('Ga: avg=%3.3f min=%3.3f max=%3.3f',mean(Ga_of_ta),min(Ga_of_ta),max(Ga_of_ta)));
    % Verifications
	verify_inhom_fix_1D(context);
	Ge  = [0,0,Ge];
    
    % SLR & windowing
    if (SLR_flag)
        declare_stage('RE-design pulse using SLR');
        sw = 200e+3;
        pulse_fname = 'D:\PhD\MATLAB\Simulations\Ver_4\chirp4SLR';
WRONG PARAMETERS PASSED TO SLR FUNCTION::  chirp_rot_SLR = convert_pulse_to_SLR(z_axis,Tp,Phi_e_sq,OmegaE,ze_of_OmegaE,chirp_rot,sw,tpwr90,pw90,RF_ESP_SLR,pulse_fname);
        chirp_rot     = chirp_rot_SLR;
        phi_chirp_rot = phase(chirp_rot_SLR);
        dte = 1/sw;
        te  = linspace(0,Tp,length(chirp_rot));
    end;
    % Windowing
	if (RFnGRAD_Win_Flag)
		Ga_of_ta  = window_grad(Ga_of_ta,ta,Ga_ESP,DEBUG_FLAG);
		if (~SLR_flag) % o/w the pulse was already smoothed
			chirp_rot = window_rf(chirp_rot,te,RF_ESP,1);%DEBUG_FLAG);
		end;
	end;
else % if (DesignPulse_Flag)
    if (SLR_flag)
		declare_substage('RE-design pulse using SLR');

% 		overlap_factor = -30/100;                                   % [%]
% 		cur_Oi = Oi - overlap_factor * dO;                          % [Hz]
% 		cur_R  = R_of_te * (1 + 2*overlap_factor);                  % [Hz/sec]
% 		cur_dO = cur_R*Tp;                                          % [Hz]
		
		for idx = 1:Nsp
			% Calculate the pulse frequency range
% 			Oi_n(idx) = cur_Oi + (idx-1)*cur_dO;                    % [Hz]
% 			Of = Oi_n(idx) + cur_dO;                                % [Hz]
% 			alpha0_n(idx) = - Tp * ((Oi_n(idx))^2) / (2*cur_dO);
% 			alpha1_n(idx) = gammaHz * Ge(3) * Tp * Of / cur_dO;
			
			Oi_n(idx) = Oi + (idx-1)*dO;                            % [Hz]
			Of = Oi_n(idx) + dO;                                    % [Hz]
			alpha0_n(idx) = - Tp * ((Oi_n(idx))^2) / (2*dO);
			alpha1_n(idx) = gammaHz * Ge(3) * Tp * Of / dO;
		end;
		alpha2 = - ((gammaHz*Ge(3))^2) / (2*R_of_te);
		
		OmegaE      = gammaHz*Ge(3)*z_axis;
		OmegaE_axis = OmegaE;
		z_of_OmegaE = interp1(OmegaE,z_axis,OmegaE);
		swHz        = 1/dte;
		isPI        = 0;
		tpwr90      = 60;
		pw90        = 23;
		if (Nsp > 1)
			n = round(linspace(0,length(z_axis),Nsp+1));
			Phi_e_sq_rad = [];
			for sp_idx = 1:Nsp
				ax = z_axis(n(sp_idx)+1 : n(sp_idx+1));                               % Spatial range of current SP
				mid_val = mean([z_axis(n(sp_idx)+1), z_axis(n(sp_idx+1))]);           % Middle value of the range
				Phi_e_sq_rad((end+1):n(sp_idx+1)) = 2*pi*alpha2*((ax - mid_val) .^ 2);  % Frequency response
			end;
		else
			Phi_e_sq_rad = 2*pi*alpha2*(z_axis.^2);
		end;

		if (DEBUG_FLAG)
		figure;
		subplot(3,1,1); plot(z_axis      ,OmegaE*1E-3 ,'.-'); title('OmegaE(z)  '); xlabel('z-axis [cm]'   ); ylabel('Omega(z) [kHz]'       ); grid;
		subplot(3,1,2); plot(OmegaE*1E-3 ,z_of_OmegaE ,'.-'); title('z(OmegaE)  '); xlabel('Omega(z) [kHz]'); ylabel('z-axis [cm]'          ); grid;
		subplot(3,1,3); plot(z_axis      ,Phi_e_sq_rad,'.-'); title('Phi e sq(z)'); xlabel('z-axis [cm]'   ); ylabel('Phi_e_sq_rad(z) [rad]'); grid;
		end;
		[chirp_rot,rfpwr,rffpwr] = convert_pulse_to_SLR(z_axis,Tp,rfwdth,Phi_e_sq_rad,OmegaE,OmegaE_axis,...
		                                                z_of_OmegaE,swHz,isPI,tpwr90,pw90,RF_ESP_SLR,DEBUG_FLAG);
		phi_chirp_rot = phase(chirp_rot);
    else % if (SLR_flag)
		[chirp_rot, phi_chirp_rot] = gen_Chirp_pulse(te,Oi,R_of_te,context);
	end;
end;

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
% Pre PI pulse purge in order to reverse acquisition direction
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
% PI pulse
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
% Post PI pulse purge in order to reverse acquisition direction
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
if (DesignPulse_Flag)
	M_final_z = abs(Sig / dZa);
else
	R_rad   = 2*pi*R_of_te;                                                    % [rad/sec^2]
	Oi_rad  = 2*pi*Oi;                                                         % [rad/sec]
	alpha_0 = (((omega_cs - Oi_rad)^2) ./ (2*R_rad)) - Tp*omega_cs;            % [rad]
	alpha_1 = gammaHz*2*pi * Ge(3) * (((omega_cs - Oi_rad) ./ R_rad) - Tp);    % [rad/cm]
	alpha_2 = ((gammaHz*2*pi)^2) * (Ge(3)^2) ./ (2*R_rad);                     % [rad/cm^2]
	
	za_axis = (gammaHz*2*pi * Ga(3) * ta - alpha_1) / (2*alpha_2);

	M_final_z = abs(Sig / sqrt(pi*j/alpha_2));
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
	plot(sample_z_axis, M_init(:,3),'k-', hi_res_za, M_final_z(end:-1:1) ,'b-');
else
	plot(sample_z_axis, M_init(:,3),'k-', za_axis, M_final_z(end:-1:1) ,'b-');
end;
title('Initial & Reconstructed Magnetization');  xlabel('z [cm]');  ylabel('Magnetization');
axis([min(sample_z_axis),max(sample_z_axis),-0.1,max(M_final_z)+0.2]); % legend({'Initial M','Reconstructed M'});
set_gca; set(gca,'FontSize',20);

% -------------------------------------------------------
%  Create output file for Multi-SP post-processing stage
% -------------------------------------------------------
if (Nsp > 1)
	fn = sprintf('D:\\PhD\\Matlab\\Simulations\\Ver_4\\UF1D90\\MultSP_sim\\%3.0f.fid',FIDnum);
% 	while (exist(fn,'file') == 2)
% 		op.Resize      = 'on';
% 		op.WindowStyle = 'normal';
% 		FIDnumCell = inputdlg('Please enter a new FID number',sprintf('File %3.0f exist',FIDnum),1,{'FID #'},op);
% 		if (isempty(FIDnumCell))
% 			FIDnum = FIDnum + 1;
% 			uiwait(msgbox(sprintf('Incrementing FIDnum by 1 to %3.0f',FIDnum)));
% 		else
% 			FIDnum = str2double(FIDnumCell{1});
% 		end;
% 		fn = sprintf('D:\\PhD\\Matlab\\Simulations\\Ver_4\\UF1D90\\MultSP_sim\\%3.0f.fid',FIDnum);
% 	end;
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
	fprintf(fd,'\n');
	fclose(fd);

% 	close all;
% 	MultSP_1DUF_pp(FIDnum);
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
