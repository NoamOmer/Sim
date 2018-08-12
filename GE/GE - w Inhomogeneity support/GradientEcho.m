% Gradient-Echo
clear all; close all; pack;
simulation_name = 'Gradient-Echo Simulation';
declare_start(simulation_name);
compile_C_code();
context = 'ge_set_globals';
set_context;

% ------------------------------------------------
% Retrieve B0 inhomogeneity and set a spatial axis
% ------------------------------------------------
[map_fn,P_DNu0,z_axis,map_tof,op,mean_std] = Inhomo_Zmap;
if (map_tof ~= tof)    warning('Map tof does not match current tof');                      end;
if (DEBUG_FLAG >= 2)   fh = figure; plot(z_axis,polyval(P_DNu0,z_axis)); hold on;          end;

z_axis = min(Zie,Zfe):dz:max(Zie,Zfe);
dz     = abs(z_axis(2)-z_axis(1));
Lz     = abs(z_axis(end)-z_axis(1));
if (DEBUG_FLAG)    figure(fh); hold on; plot(z_axis,polyval(P_DNu0,z_axis),'r'); grid on;  end;
disp(sprintf('New z-axis range:\n dz=%d\n Lz=%d',dz,Lz));

if (DEBUG_FLAG >= 2)
    figure(fh); hold on; plot(z_axis,polyval(P_DNu0,z_axis),'c--')
    title('Field Inhomogeneity'); xlabel('z-axis [cm]'); ylabel('Inhomogeneity [Hz]');
    legend({'Initial','After range set','After Inhomo scaling'},'Location','Best'); grid;
end;

% Now that the z_axis has setteled -- define the inhomogeneity
dOmega_0 = polyval(P_DNu0,z_axis);
dB0z = (dOmega_0/gammaHz)*1E-4;

% -------------------------------------------------------------------------
% Create a shpaed sample
% -------------------------------------------------------------------------
sample = ones(1,length(z_axis));
M_init = transpose(sample) * M0;

if (DEBUG_FLAG)    figure;
plot(z_axis,M_init(:,3),'.-'); title('Sample(z)'); xlabel('z-axis [cm]'); ylabel('Magnetization [amper/m]');
set_gca;           end;

% -------------------------------------------------------------------------
% Excitation (non slice-selective hard pulse)
% -------------------------------------------------------------------------
tic
declare_stage('Excitation');
[B_RF_rot, B_RF_t] = gen_RF_pulse(Rect_P, flip_a*pi/180, 0, omega_cs, context);
disp(sprintf('%d Excitation duration = %5.2f [us]',flip_a,(B_RF_t(end))*1E+6));

B_eff_rot = B_RF_rot + (transpose(ones(1,length(B_RF_t)))) * B_rot;
if (C_CODE)
 	Exc_Gradient = zeros(1,length(B_eff_rot(:,1)));
  	st = set_evolve_M_n_acq1Dgrad_t_CPP_struct(M_init,B_eff_rot,dt,dB0z,Exc_Gradient,z_axis,...
  	                                           RH_flag,0,Inhomo_Flag,Relax_Flag,T1,T2,context);
 	[Mx,My,Mz,dummy1,dummy2] = evolve_M_n_acq1Dgrad_t_CPP(st);
 	M = [transpose(Mx),transpose(My),transpose(Mz)];
else
    Exc_Gradient = [0,0,0];
    M = evolve_M_n_acq(M_init,B_eff_rot,dB0z,Exc_Gradient,dt,z_axis,0,Inhomo_Flag,0);
end;
plot_M_n_Phase(0,M_init,M,'M(z) after pi/2 excitation pulse','Phase after pi/2 excitation',{'\phi(z)'},context);
toc

% -------------------------------------------------------------------------
% Pre-acquisition refocusing gradient
% -------------------------------------------------------------------------
declare_stage('Pre-acq refocusing gradient'); tic
B_eff_rot = B_rot;
if (C_CODE)
	st = set_evolve_M_CPP_struct(M,B_eff_rot,Tpr,dB0z,Gpr,z_axis,RH_flag,...
	                             Inhomo_Flag,Relax_Flag,T1,T2,context);
	[Mx,My,Mz] = evolve_M_CPP(st);
 	M = [transpose(Mx),transpose(My),transpose(Mz)];
else
    M = evolve_M(M,B_eff_rot,dB0z,[0,0,Gpr],Tpr,z_axis,Inhomo_Flag);
end;
plot_M_n_Phase(0,M_init,M,'M(z) after Pre-acq refocusing gradient','Phase after Pre-acq refocusing gradient',{'\phi(z)'},context);
toc

% -------------------------------------------------------------------------
% Acquisition
% -------------------------------------------------------------------------
declare_stage('Acquisition'); tic
B_eff_rot = transpose(ones(1,length(ta))) * B_rot;
if (C_CODE)
    Exc_Gradient = Ga * ones(1,length(ta));
    st = set_evolve_M_n_acq1Dgrad_t_CPP_struct(M,B_eff_rot,dta,dB0z,Exc_Gradient,z_axis,RH_flag,1,...
                                               Inhomo_Flag,Relax_Flag,T1,T2,context);
    [Mx,My,Mz,sig_real,sig_imag] = evolve_M_n_acq1Dgrad_t_CPP(st);
    M = [transpose(Mx),transpose(My),transpose(Mz)];
    Sig = transpose(sig_real) + i*transpose(sig_imag);
else
    Exc_Gradient = [0,0,1]*Ga;
    [M,Sig] = evolve_M_n_acq(M,B_eff_rot,dB0z,Exc_Gradient,dta,z_axis,1,Inhomo_Flag,0);
    sig_real = real(Sig);
    sig_imag = imag(Sig);
end;

if (DEBUG_FLAG) figure;
plot(ta*1E+3,sig_real,'b-',ta*1E+3,sig_imag,'m-',ta*1E+3,abs(Sig),'k-');
title('FID'); xlabel('Time [ms]'); ylabel('Signal'); legend({'Real','Image','Abs'});
set_gca; end;
toc

% -------------------------------------------------------------------------
% Post processing - FFT & display
% -------------------------------------------------------------------------
declare_stage('Post-Processing');
%sig_pad = zero_pad(transpose(Sig));
%dta = (ta(2)-ta(1));
%ta = ta(1) : dta : dta*length(sig_pad);

fft_sig = fftshift(fft(Sig));
figure;
if (RH_flag)
    plot(nu_axis(end:-1:1),-imag(fft_sig),'k.-');  % Imaginary part gets a minus sign when using RHR
else
    plot(nu_axis          ,+imag(fft_sig),'k.-');
end;    
title(sprintf('Sig FT'));
xlabel('Frequency [Hz]'); ylabel('FT(Signal)');  set_gca;

disp(sprintf(['\nat            = %9.1f [ms]  '  ...
              '\nsw            = %9.3f [kHz] '  ...
              '\nnp            = %9.1f       '  ...
              '\ntof           = %9.3f [Hz]  '  ...
              '\nflip-angle    = %9.1f [deg] '  ...
              '\nZie           = %9.3f [cm]  '  ...
              '\nZfe           = %9.3f [cm]  '  ...
              '\nLz            = %9.3f [cm]  '  ...
              '\ndz            = %9.5f [mm]  '  ...
              '\nRelaxation    = %9.1f       '  ...
              '\nT1            = %9.1f [ms]  '  ...
              '\nT2            = %9.1f [ms]  '  ...
              '\nRH_Flag       = %9.0f       '  ...
              '\nmap_fn        = %s          '  ...
              '\nmap_op        = %9.0f       '],...
              Ta*1e+3,sw*1e-3,n_pts,tof,flip_a,Zie,Zfe,Lz,dz*10,...
              Relax_Flag,T1*1E+3,T2*1E+3,RH_flag,map_fn,op));

declare_end(simulation_name);
return;

