% 1 Dimensional Ultrafast MRI using 180-chirp pulse
% 1D UF Version 1
% 1. Excitation - 90 hard pulse - non slice-selective - inhomogeneity: off
% 2. PI chirp - inhomogeneity: on
% 3. Purge - inhomogeneity: off
% 4. Acquisition

% 1D UF Version 2
% Balance the inhomogeneity evolution during this stage by the selective excitation purge stage which
% occurs prior the 180 pulse, thus being of the same duration gains the same phase, yet with opposite sign.
% 1. Excitation - 90 hard pulse (T1) - non slice-selective - inhomogeneity: on
% 2. Free evolution - inhomogeneity: on (T2)
% 3. PI chirp - inhomogeneity: on
% 4. Purge - inhomogeneity: on - duration: 0.5*T1 + T2 - balancing the pre-chirp inhomogeneity evolution
% 5. Acquisition
clear all; %close all; pack;
version = 2;
mex 'D:\PhD\MATLAB\Simulations\Ver_4\C code\evolve_M_n_acq1Dgrad_t_CPP.c'
mex 'D:\PhD\MATLAB\Simulations\Ver_4\C code\evolve_M_CPP.c'
context = 'UF_1D180_set_globals';
set_context;
declare_start(sim_name);  tic

% -------------------------------------------------------------------------
% Retrieve the B0 inhomogeneity and set the spatial axis
% -------------------------------------------------------------------------
declare_stage('Retrieve the B0 inhomogeneity and set the spatial axis');
% [map_fn,P_DNu0,z_axis,tof,op,mean_std] = Inhomo_Zmap;
if (DEBUG_FLAG >= 1)
    figure;
	plot(z_axis,polyval(P_DNu0,z_axis),'.-'); hold on;
    title('Field inhomogeneity'); xlabel('z-axis [cm]'); ylabel('\Delta\Nu [Hz]');
end;

% Inhomogeneous field tests
P_DNu0 = Map_Mult_Factor*P_DNu0;
if (DEBUG_FLAG >= 1)    plot(z_axis,polyval(P_DNu0,z_axis),'c');    end;

% Interpolate the internal z_axis so that the RF dwell time in the te axis
% calculation will match the global dwell-time value
n_RF_pts  = 1 + Tp/RFdt;
z_axis = linspace(Zie,Zfe,n_RF_pts*Ze_INTERP_FACTOR);
dz = z_axis(2) - z_axis(1);
Lz = abs(z_axis(end) - z_axis(1));
disp(sprintf('UF1D180:\n dz=%d;\n Lz=%d',dz,Lz));
if (DEBUG_FLAG >= 1)
    plot(z_axis,polyval(P_DNu0,z_axis),'m.');
    legend({'Initial','After Inhomo scaling','After z-interp and range setting'},'Location','Best'); grid;
end;

% -------------------------------------------------------------------------
% Create a shaped sample
% -------------------------------------------------------------------------
% Create a new z-axis for the the sample
sample_z_axis = (sample_Zi : (z_axis(2) - z_axis(1)) : sample_Zf);

% sample = exp(-2E+1*((z_axis-00*dz).^2)) + exp(-2E+1*((z_axis-40*dz).^2)); % Smooth shape (gaussian) decreases the effect of Gibbs phenomenon
% sample = exp(-1E+5*((z_axis*dz).^2)) + exp(-1E+5*((z_axis*dz).^2));
% % Square
% sample = [zeros(1,750), ones(1,length(sample_z_axis)-1500), zeros(1,750)];   % This shape is convenient for examining the excitation efficiency
% % Triangle
% F = 500;
% sample = [zeros(1,F), 1:(length(z_axis)/2-F), (length(z_axis)/2-F):-1:1, zeros(1,F)]; sample = sample/max(sample);
% if (length(sample) < length(z_axis))   sample(end+1) = sample(end);   end;
% Ones
sample = ones(1,length(sample_z_axis));
M_init = transpose(sample) * M0;
if (DEBUG_FLAG >= 2)
	figure; plot(sample_z_axis,M_init(:,3),'.-'); title('Sample(z)'); xlabel('z-axis [cm]'); ylabel('M'); set_gca;
end;

% Now that the z_axis has setteled -- define the inhomogeneity
B_lab = B0_lab + [0, 0, (omega_cs /gamma_T)];
B_rot = B_lab  - [0, 0, (omega_rot/gamma_T)];
dOmega_0 = polyval(P_DNu0,sample_z_axis);     % [Hz]
dB0z = (dOmega_0/gammaHz)*1E-4;               % [T]

% -------------------------------------------------------------------------
%  Excitation (non slice-selective hard pulse)
% -------------------------------------------------------------------------
declare_stage('Excitation (non slice--selective hard pulse)');
[B_RF_rot, B_RF_t] = gen_RF_pulse(Rect_P, pi/2, 0, omega_cs, context);
disp(sprintf('pi/2 Excitation duration = %d [us]',(B_RF_t(end))*1E+6));
B_eff_rot = B_RF_rot + (transpose(ones(1,length(B_RF_t)))) * B_rot;

switch (version)
case 1
    disp('**  Version 1: simple hard pulse  **');
    acq_flag = 0;
    [M, dummy]=evolve_M_n_acq(M_init,B_eff_rot,zeros(1,length(dB0z)),[0,0,0],dt,sample_z_axis,...
                              acq_flag,Inhomo_Flag,DEBUG_FLAG);
    plot_M_n_Phase(0,M_init,M,'M(z) after pi/2 excitation pulse','Phase after pi/2 excitation',{'\phi(z)'},context);
    
case 2
    disp('**  Version 2: hard pulse + delay **');
    % Hard pulse
    if (C_CODE)
        B_eff_rot = [B_eff_rot(1,1) B_eff_rot(1,2) B_eff_rot(1,3)];
        zGrad = 0;
        st = set_evolve_M_CPP_struct(M_init,B_eff_rot,B_RF_t(end),dB0z,zGrad,sample_z_axis,RH_flag,...
                                     Inhomo_Flag,Relax_Flag,T1,T2,context);
        [Mx,My,Mz] = evolve_M_CPP(st);
        M = [transpose(Mx),transpose(My),transpose(Mz)];
    else
        acq_flag = 0;
        zGrad    = [0,0,0];
        [M, dummy]=evolve_M_n_acq(M_init,B_eff_rot,dB0z,zGrad,dt,sample_z_axis,acq_flag,Inhomo_Flag,DEBUG_FLAG);
    end;
    plot_M_n_Phase(0,M_init,M,'M(z) after pi/2 excitation pulse',...
        'Phase after pi/2 excitation (should be const. as long as inhomogeneity --> 0',...
        {'\phi(z)'},context);

    % Delay
	B_eff_rot = B_rot;
	Exc_Delay = Tpr - B_RF_t(end);
	disp(sprintf('Post-excitation delay = %d [us]',Exc_Delay*1E+6));
	if (C_CODE)
		st = set_evolve_M_CPP_struct(M,B_eff_rot,Exc_Delay,dB0z,0,sample_z_axis,RH_flag,...
		                             Inhomo_Flag,Relax_Flag,T1,T2,context);
		[Mx,My,Mz] = evolve_M_CPP(st);
		M = [transpose(Mx),transpose(My),transpose(Mz)];
	else
		M = evolve_M(M,B_eff_rot,dB0z,[0,0,0],Exc_Delay,sample_z_axis,Inhomo_Flag);
	end;
	plot_M_n_Phase(0,M_init,M,'M(z) after post-excitation delay','Phase after post-excitation delay',{'\phi(z)'},context);

end;

% -------------------------------------------------------------------------
% Pre 180 purge
% -------------------------------------------------------------------------
declare_stage('Pre PI-Chirp purge');
B_eff_rot = B_rot;
Pre_180_Exc_Delay = 2.5e-3;  % [sec]
Pre_180_zGrad     = 40;      % [G/cm]
disp(sprintf('Pre 180 purge: %d [us] %d [G/cm]',Pre_180_Exc_Delay*1E+6,Pre_180_zGrad));
st = set_evolve_M_CPP_struct(M,B_eff_rot,Pre_180_Exc_Delay,dB0z,Pre_180_zGrad,sample_z_axis,RH_flag,...
                             Inhomo_Flag,Relax_Flag,T1,T2,context);
[Mx,My,Mz] = evolve_M_CPP(st);
M = [transpose(Mx),transpose(My),transpose(Mz)];

% -------------------------------------------------------------------------
% Design & Apply PI-Chirp pulse
% -------------------------------------------------------------------------
declare_stage('Design PI-Chirp pulse');
if (Inhomo_Fix_Flag == 0)
    P_DNu0_ = P_DNu0;
    P_DNu0  = P_DNu0*0;
end;

if (SLR)
	OmegaE       = gammaHz*UD.GPEe*z_axis_180;
	OmegaE_axis  = OmegaE;
	z_of_OmegaE  = interp1(OmegaE,z_axis_180,OmegaE);
	Phi_e_sq_rad = 2*pi*alpha2_180*(z_axis_180.^2);
	swHz         = 1/dtse;
	isPI         = 1;
	if (DEBUG_FLAG >= 2)
	figure;
	subplot(3,1,1); plot(z_axis_180 ,OmegaE*1E-3 ,'.-'); title('OmegaE(z)  '); xlabel('z-axis [cm]'   ); ylabel('Omega(z) [kHz]'       ); grid;
	subplot(3,1,2); plot(OmegaE*1E-3,z_of_OmegaE ,'.-'); title('z(OmegaE)  '); xlabel('Omega(z) [kHz]'); ylabel('z-axis [cm]'          ); grid;
	subplot(3,1,3); plot(z_axis_180 ,Phi_e_sq_rad,'.-'); title('Phi e sq(z)'); xlabel('z-axis [cm]'   ); ylabel('Phi_e_sq_rad(z) [rad]'); grid;
	end;
	disp(sprintf('Designing SE SLR pulse. Tse=%3.3f, rf180wdth=%3.3f',UD.Tse,UD.rf180wdth));
	[chirp180,rf180pwr,rf180fpwr] = convert_pulse_to_SLR(z_axis_180,UD.Tse,UD.rf180wdth,Phi_e_sq_rad,...
	                                                     OmegaE,OmegaE_axis,z_of_OmegaE,swHz,isPI,...
	                                                     UD.tpwr90,UD.pw90,UD.RF180_ESP,0);
else
	if RH_flag
		[te,dte,O_of_te,Ge,Gpr,Ga_of_ta,R_of_te,chirp_rot,phi_chirp_rot] = design_inhom_fix_UF_1D180_RHR(context);
	else
		[te,dte,O_of_te,Ge,Gpr,Ga_of_ta,R_of_te,chirp_rot,phi_chirp_rot] = design_inhom_fix_UF_1D180(context);
	end;
	if (Inhomo_Fix_Flag == 0)
		P_DNu0 = P_DNu0_;
	end;
	DOkHz = 1E-3*abs(O_of_te(end) - O_of_te(1));
	disp(sprintf('Ga: avg=%3.3f min=%3.3f max=%3.3f',mean(Ga_of_ta),min(Ga_of_ta),max(Ga_of_ta)));
	% Verifications
	verify_inhom_fix_1D180(context);
end;

% Windowing
if (RFnGRAD_Win_Flag)
	Ga_of_ta  = window_grad(Ga_of_ta,ta,Ga_ESP,DEBUG_FLAG);
	chirp_rot = window_rf(chirp_rot,te,RF_ESP,DEBUG_FLAG);
end;

Ge3D  = [0,0,Ge];
Gpr3D = [0,0,Gpr];

x_chirp_rot = real(chirp_rot);
y_chirp_rot = imag(chirp_rot);
z_chirp_rot = zeros(1,length(chirp_rot));

if (DEBUG_FLAG >= 2)
	figure;   plot(te*1e+3, abs(chirp_rot)*gamma_T/2/pi, 'k.-');
	title('Chirp pulse amplitude in rotating frame');  xlabel('Time [ms]');  ylabel('Apmlitude [Hz]'); set_gca;

	figure;   plot(te*1e+3, x_chirp_rot, 'k.-', te*1e+3, y_chirp_rot, 'b.-');
	title('Chirp pulse in rotating frame');  xlabel('Time [ms]');  ylabel('Magnetic field'); 
	legend({'x-coordinate','y-coordinate'});   set_gca;

	figure;   plot(te*1e+3, phi_chirp_rot, 'b.-');
	title('Chirp angle in rotating frame');  xlabel('Time [ms]');  ylabel('Angle [rad]'); set_gca;
end;

declare_stage('Apply PI-Chirp pulse');
B_eff_rot = [transpose(x_chirp_rot),transpose(y_chirp_rot),transpose(z_chirp_rot)] + (transpose(ones(1,length(chirp_rot))))*B_rot;
if (C_CODE)
	zGrad = ones(1,length(B_eff_rot(:,1)))*Ge3D(3);
	st = set_evolve_M_n_acq1Dgrad_t_CPP_struct(M,B_eff_rot,dte,dB0z,zGrad,sample_z_axis,RH_flag,...
	                                           0,Inhomo_Flag,Relax_Flag,T1,T2,context);
	[Mx,My,Mz,dummy1,dummy2] = evolve_M_n_acq1Dgrad_t_CPP(st);
	M = [transpose(Mx),transpose(My),transpose(Mz)];
else
    [M, dummy] = evolve_M_n_acq(M, B_eff_rot, dB0z, Ge3D, dte, sample_z_axis, 0, Inhomo_Flag, 0);
end;
plot_M_n_Phase(0,M_init,M,'M(z) after RF Chirp pulse','Phase after RF Chirp pulse',{'\phi(z)'},context);

% -------------------------------------------------------------------------
% Post 180 purge
% -------------------------------------------------------------------------
declare_stage('Post PI-Chirp purge');
B_eff_rot = B_rot;
Exc_Delay = Pre_180_Exc_Delay;  % [sec]
zGrad     = Pre_180_zGrad;      % [G/cm]
disp(sprintf('Post 180 purge: %d [us] %d [G/cm]',Exc_Delay*1E+6,zGrad));
st = set_evolve_M_CPP_struct(M,B_eff_rot,Exc_Delay,dB0z,zGrad,sample_z_axis,RH_flag,Inhomo_Flag,Relax_Flag,T1,T2,context);
[Mx,My,Mz] = evolve_M_CPP(st);
M = [transpose(Mx),transpose(My),transpose(Mz)];

% -------------------------------------------------------------------------
% Purge stage
% -------------------------------------------------------------------------
declare_stage('Purge stage');
B_eff_rot = B_rot;
if (version == 1)
	disp('**  Version 1: Purge w/o inhomogeneity evolution  **');
	dB0z_purge = zeros(1,length(dB0z));
elseif (version == 2)
	disp('**  Version 2: Purge w/  inhomogeneity evolution  **');
	dB0z_purge = dB0z;
end;
if (C_CODE)
	st = set_evolve_M_CPP_struct(M,B_eff_rot,Tpr,dB0z_purge,Gpr3D(3),sample_z_axis,RH_flag,...
	                             Inhomo_Flag,Relax_Flag,T1,T2,context);
	[Mx,My,Mz] = evolve_M_CPP(st);
	M = [transpose(Mx),transpose(My),transpose(Mz)];
else
	M = evolve_M(M,B_eff_rot,dB0z_purge,Gpr3D,Tpr,sample_z_axis,Inhomo_Flag);
end;
h1 = plot_M_n_Phase(0,M_init,M,'M(z) after purge stage','\phi(z) after purge stage',{'\phi(z)'},context);

% -------------------------------------------------------------------------
% Acquire signal
% -------------------------------------------------------------------------
declare_stage('Acquiring signal');
% TEST >>>>
if (0)
	for idx=1:length(sample_z_axis)
		r = rot_rhr([0,0,1],theoretical_phi_pre_acq1(idx));
		M1(idx,1:3) = M_init(idx,3) .* transpose(r*transpose([1,0,0]));

		r2 = rot_rhr([0,0,1],theoretical_phi_pre_acq2(idx));
		M2(idx,1:3) = M_init(idx,3) .* transpose(r2*transpose([1,0,0]));

		r3 = rot_rhr([0,0,1],theoretical_phi_pre_acq3(idx));
		M3(idx,1:3) = M_init(idx,3) .* transpose(r3*transpose([1,0,0]));
	end;

	plot_M_n_Phase(h1,M_init,M,'M(z) after artificial phase setting',...
	                           '\phi(z) Before and after artificial phase setting',...
	                           {'Before','After'},...
	                           context);
end;
% <<<< End Test

B_eff_rot = transpose(ones(1,length(ta))) * B_rot;
if (C_CODE)
	st = set_evolve_M_n_acq1Dgrad_t_CPP_struct(M,B_eff_rot,dta,dB0z,Ga_of_ta,sample_z_axis,RH_flag,...
	                                           1,Inhomo_Flag,Relax_Flag,T1,T2,context);
	[Mx,My,Mz,sig_real,sig_imag] = evolve_M_n_acq1Dgrad_t_CPP(st);
	M = [transpose(Mx),transpose(My),transpose(Mz)];
	Sig = transpose(sig_real) + i*transpose(sig_imag);
else
	[M, Sig] = evolve_M_n_acq1Dgrad_t(M, B_eff_rot, dB0z, Ga_of_ta, dta, sample_z_axis, 1, Inhomo_Flag, 0);
end;

% -------------------------------------------------------------------------
% Reconstruct signal
% -------------------------------------------------------------------------
M_final_z = abs(Sig / dZa);

figure; hold;
plot(sample_z_axis, M_init(:,3),'k-', hi_res_za(end:-1:1), M_final_z ,'b-');
title('Initial & Reconstructed Magnetization');  xlabel('z [cm]');  ylabel('Magnetization');
legend({'initial','reconstructed'});
declare_end(sim_name); toc

return;

% ---------
%   JUNK
% ---------
%1% Replaced the following with a simpler and shorter code
% x_B_RF_rot = B_RF_rot(:,1);
% y_B_RF_rot = B_RF_rot(:,2);
% z_B_RF_rot = B_RF_rot(:,3);
% B_eff_rot = [x_B_RF_rot, y_B_RF_rot, z_B_RF_rot] + (ones(1,length(B_RF_t))') * B_rot;

% % Extrapolate the excitation range -- for simulation purposes
% EXTRAP_FACTOR = 1.7;
% dz_tmp = z_axis(2) - z_axis(1);
% z_axis = ((z_axis(1)*EXTRAP_FACTOR):dz_tmp:(z_axis(end)*EXTRAP_FACTOR));
% 
% % Set the excitation range (we might want to excite only a subset of the inhomogeneity range)
% Zie_idx = max(find((z_axis - Zie) < 0));
% Zfe_idx = min(find((z_axis - Zfe) > 0));
% if (isempty(Zie_idx))   Zie_idx = 1;               end;
% if (isempty(Zfe_idx))   Zfe_idx = length(z_axis);  end;
% z_axis = z_axis(Zie_idx:sign(Zfe_idx-Zie_idx):Zfe_idx);
% plot(z_axis,polyval(P_DNu0,z_axis),'r.-');
% 
% % Inhomogeneous field tests
% P_DNu0 = Map_Mult_Factor*P_DNu0;
% plot(z_axis,polyval(P_DNu0,z_axis),'c--');
% 
% % Interpolate the internal z_axis so that the RF dwell time in the te axis
% % calculation will match the global dwell-time value
% n_RF_pts  = 1 + Tp/RFdt;
% z_axis = linspace(z_axis(1),z_axis(end),n_RF_pts*2.9);
% dz = z_axis(2) - z_axis(1);
% Lz = abs(z_axis(end) - z_axis(1));
% if (DEBUG_FLAG)  disp(sprintf('UF1D180:\n dz=%d;\n Lz=%d',dz,Lz));  end;
% plot(z_axis,polyval(P_DNu0,z_axis),'m.-');
% legend({'Initial','After range set','After Inhomo scaling','After z-interp'},'Location','Best'); grid;
% 
