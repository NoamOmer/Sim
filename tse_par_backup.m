% ___________________________________________________________________________
% ===========================================================================
% SE_MC Simulation
% Rev 1.0  |  Aug 2012  |  Noam Ben-Eliezer  |  Initial version
% Rev 2.0  |  Oct 2012  |  Noam Ben-Eliezer  |  Add support for multi-T2
% Rev 2.1  |  Dec 2012  |  Noam Ben-Eliezer  |  Add support for multi-B1
% ___________________________________________________________________________
% ===========================================================================

% clcl;
clear all; %close all;
simulation_name = 'TSE Simulation';
declare_start(simulation_name);
compile_C_code();

% % % % 		global required_res;
% % % % 		reqres_idx = 0;
% % % % 		fh1=0;
% % % % % 		lgd1 = {'0.02'};
% % % % 		lgd1 = {'0.005','0.01','0.02','0.03','0.04'};
% % % % % 		for required_res = [0.02];
% % % % 		for required_res = [0.005 0.01 0.02 0.03 0.04];
% % % % 		reqres_idx = reqres_idx + 1;

% % % % % % % % % % 	global refoc_angle;
% % % % % % % % % % 	refoc_angle_idx = 0;
% % % % % % % % % % 	fh1=0;
% % % % % % % % % % 	lgd1 = {'150','180','210'};
% % % % % % % % % % 	for refoc_angle = [150,180,210]
% % % % % % % % % % 	refoc_angle_idx = refoc_angle_idx + 1;


% context = 'tse_set_globals_2013_06_10_bw352_12_12_96';
% context = 'tse_set_globals_2012_10_09';
% context = 'tse_set_globals_2012_12_26_TF10';
% context = 'tse_set_globals_2013_02_05_MID40ii';                  % MnCl2 phantom
% context = 'tse_set_globals_2013_02_05_MID40ii_forNoiseTest';     % MnCl2 phantom
% context = 'tse_set_globals_2013_03_19_TF6_paperRev';
% context = 'tse_set_globals_2013_03_19_T1';
% context = 'tse_set_globals_2013_07_24_TE12_RFnorm';
% context = 'tse_set_globals_2013_07_30_MID472_3_7_8';
% context = 'tse_set_globals_2013_07_30_MID479';
% context = 'tse_set_globals_2013_09_15';
% context = 'tse_set_globals_2013_08_07';                          % Joey Prostate (SEMC comparison to RAISE + CS)
% context = 'tse_set_globals_2014_02_03';                          % Giuseppe C. Hip Tim Trio
% context = 'tse_set_globals_2014_Apr_oldSEMCdata';
% context = 'tse_set_globals_2014_05_18_MID153s';
% context = 'tse_set_globals_2014_06_09';
% context = 'tse_set_globals_Hip_surgical_Skyra_TE12';
% context = 'tse_set_globals_2014_London';
% context = 'tse_set_globals_Akio';
% context = 'tse_set_globals_2014_10_Ileana';
context = 'tse_set_globals';


tic;
set_context;
% profile on;


% =========================================================================
% Create a shaped sample
% =========================================================================
[sample_z_axis,M_init] = UD1D_create_sample(context);

%
%  === Initialize refocusing pulse ===
%
if load_ext_RF
    [b1_refoc_SiemensPTA, N_refoc] = read_siemens_PTA_pulse(Refname);
    dte_refoc = Trefoc/N_refoc;
else
	b1_refoc_SiemensPTA  = [];
    dte_refoc = Trefoc/N_refoc;
    pulse_prop.type = Sinc_P;
    pulse_prop.n_lobes = n_lobes_refoc;
    pulse_prop.Tp = Trefoc;
    pulse_prop.calib = 0;
    pulse_prop.Hanning_Win_Flag = Henning_Flag;              % Because that's what Siemens does to it's SINC pulse...
	Refname = sprintf('%1.0 lobes SINC',pulse_prop.n_lobes);
end;

% =========================================================================
% Choose pre-phasing scheme -- note its effect on the TE timing calculation
% =========================================================================
[Gpr0,Gpr1,Gpr2,Grew] = tse_calc_prephase_scheme(SE_flag,pre_phas_scheme,Gpr);


% =========================================================================
% TE timing Preparation -- Follow Siemens pulse sequence scheme
% =========================================================================
[TE_fill_timePre,TE_fill_timePst,Tpr0] = tse_prep_TE_timing(TE_arr,Texc,TexcRefocus,Tpr,Gpr0,Grew,Trefoc,Tdly,Ta,ETL);


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
B0_N = length(B0_scaling_arr);
B1_N = length(B1_scaling_arr);
T1_N = length(T1_tse_arr);
T2_N = length(T2_tse_arr);
% echo_train_modulation = zeros(B0_N,B1_N,T1_N,T2_N,ETL);
B0_idx = 1:B0_N;
B1_idx = 1:B1_N;
T1_idx = 1:T1_N;
T2_idx = 1:T2_N;
%             1       2       3       4
I = allcomb(B0_idx, B1_idx, T1_idx, T2_idx);
I1 = I(:,1);
I2 = I(:,2);
I3 = I(:,3);
I4 = I(:,4);
N_idx = B0_N*B1_N*T1_N*T2_N;
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%                                                       Start  T1, T2 and B1 Loops
echo_train_modulation = zeros(N_idx,ETL);

% [T2_tse_arr,T1_tse_arr,exc_angle,refoc_angle] = init_parfor(T2_tse_arr,T1_tse_arr,exc_angle,refoc_angle);

parfor idx_ems = 1:N_idx,
    T2_tse          =                   T2_tse_arr(I4(idx_ems));
    T1_tse          =                   T1_tse_arr(I3(idx_ems));
    exc_angle_tse   = exc_angle   * B1_scaling_arr(I2(idx_ems));
    refoc_angle_tse = refoc_angle * B1_scaling_arr(I2(idx_ems));
    
    dOmega_0 = B0_scaling_arr(I1(idx_ems))*ones(1,length(sample_z_axis)); % [Hz]
    dB0z     = (dOmega_0/gammaHz)*1E-4;                                   % [T]
    % dB0z     = zeros(1,length(sample_z_axis));
    
%     fprintf('TSE_fn = %s\n\n',TSE_fn);
    
    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
% -------------------------------------------------------------------------
    % RF Pulses preparation
% -------------------------------------------------------------------------
%  === prepare excitation pulse ===
	[dte_exc,b1_Exc,B_Exc_t,GverseExc] = tse_prep_exc_RF(Excfn,Texc,exc_angle_tse,exc_phase,isVERSEexc,...
	                                                     TexcVERSE1,TexcVERSE2,GexcVERSE1,GexcVERSE2);

												 
												 
%  === prepare refocusing pulse ===
	[B_refoc,B_refoc_minus,B_refoc_180,GverseRef] = tse_prep_refoc_RF(load_ext_RF,b1_refoc_SiemensPTA,pulse_prop,refoc_angle_tse,refoc_phase,Trefoc,dte_refoc,...
	                                                                  isVERSEref,TrefVERSE1,TrefVERSE2,GrefVERSE1,GrefVERSE2);
							  
	
    % =========================================================================
    % =========================================================================
    %  Excitation
    % =========================================================================
    % =========================================================================
%     declare_stage('Excitation');
    
    % --------------------------------------------------
    %  Apply excitation pulse
    % --------------------------------------------------
    % b1_Exc = b1_Exc*0;
    % B_RF_rot(:,1) = real(b1_Exc);
    % B_RF_rot(:,2) = imag(b1_Exc);
    %     B_RF_rot(:,3) = zeros(length(b1_Exc),1);
    B_RF_rot = cat(2,real(b1_Exc(:)), imag(b1_Exc(:)), zeros(length(b1_Exc),1));
    B_eff_rot = B_RF_rot + (transpose(ones(1,length(B_Exc_t)))) * B_rot;
    
    if (isVERSEexc)
        Exc_Gradient = GverseExc;
    else
        Exc_Gradient = Gexc*ones(1,length(B_eff_rot(:,1)));
    end;
    
    st = set_evolve_M_n_acq1Dgrad_t_CPP_struct(M_init,B_eff_rot,dte_exc,dB0z,Exc_Gradient,sample_z_axis,...
        RH_flag,0,Inhomo_Flag,Relax_Flag,T1_tse,T2_tse,context);
    [Mx,My,Mz,dummy1,dummy2] = evolve_M_n_acq1Dgrad_t_CPP(st);
    M = [transpose(Mx),transpose(My),transpose(Mz)];
    plot_M_n_Phase(0,sample_z_axis,M_init,M,'M(z) after pi/2 excitation pulse','Phase after pi/2 excitation',{'\phi(z)'});
    % grads{1} = GradPulse;
    % grads{1}.CalcAmpl(StartTime,TotalTime,Amplitude,Axis);   % StartTime [us];
    
   
    % --------------------------------------------------
    %  1. Refocus excitation gradient
    %  2. Apply first crusher
    %  3. Apply pre TSE loop RO pre-phasing gradient
    %  (assumes excitation pulse was not an inversion pulse)
    % --------------------------------------------------
    if (isVERSEexc)
        GexcRefocus = -0.5*(sum(GverseExc*dte_exc)) / TexcRefocus;
    else
        GexcRefocus = -0.5*(Gexc*Texc)              / TexcRefocus;
    end;
    GcrushLocal = Gcrush*Tcrush                 / TexcRefocus;
    Gpr0Local   = Gpr0*Tpr0                     / TexcRefocus;
    
    % Verify the excitation refocusing gradient value
    if (GexcRefocus_exp ~= 0)
        if (abs(abs(GexcRefocus_exp) - abs(GexcRefocus+GcrushLocal)) > 1e-1)
            fprintf('WARNING: Simualted excitation refocusing gradient (%3.3f) does not match exp. value (%3.3f)\n',(GexcRefocus+GcrushLocal),GexcRefocus_exp);
        else
            fprintf('Simualted excitation refocusing gradient (%3.3f) matches exp. value (%3.3f)\n',(GexcRefocus+GcrushLocal),GexcRefocus_exp);
        end;
    end;
    
    B_eff_rot = B_rot;
    st = set_evolve_M_CPP_struct(M,B_eff_rot,TexcRefocus,dB0z,(GexcRefocus+GcrushLocal+Gpr0Local),sample_z_axis,RH_flag,...
        Inhomo_Flag,Relax_Flag,T1_tse,T2_tse,context);
    [Mx,My,Mz] = evolve_M_CPP(st);
    M = [transpose(Mx),transpose(My),transpose(Mz)];
    plot_M_n_Phase(0,sample_z_axis,M_init,M,'M(z) after excitation refocusing gradient','Phase after excitation refocusing gradient',{'\phi(z)'});
    
   
    % =========================================================================
    % TSE loop
    % =========================================================================
%     declare_stage('TSE loop');
    % Sig      = [];
    % Im       = [];
    % Mxy_all  = [];
    % Mz_all   = [];
    
    % % --------------------------------------------------
    % % pre  TSE loop RO pre-phasing gradient (moved to be part of post-excitation gradient event)
    % % --------------------------------------------------
    % if (Gpr0 ~= 0)
    % 	B_eff_rot = B_rot;
    % 	st = set_evolve_M_CPP_struct(M,B_eff_rot,Tpr0,dB0z,Gpr0,sample_z_axis,RH_flag,Inhomo_Flag,Relax_Flag,T1_tse,T2_tse,context);
    % 	[Mx,My,Mz] = evolve_M_CPP(st);
    % 	M = [transpose(Mx),transpose(My),transpose(Mz)];
    % end;
    
    % --------------------------------------------------
    % TSE loop
    % --------------------------------------------------
    integrated_loc = [];
    echo_amp       = zeros(1,ETL);
    fprintf('B0 idx = %3.0d(%3.0d)    B1 idx = %3.0d(%3.0d)    T1 idx = %3.0d(%3.0d)    T2 idx = %3.0d(%3.0d)   idx_ems = %3.0d(%3.0d)\n',...
        I1(idx_ems),B0_N,I2(idx_ems),B1_N,I3(idx_ems),T1_N,I4(idx_ems),T2_N,idx_ems,N_idx);
    for tse_idx = 1:ETL
        % 	fprintf('B1 idx = %3.0f(%3.0f)    T1 idx = %3.0f(%3.0f)    T2 idx = %3.0f(%3.0f)    tse idx = %3.0f(%3.0f)\n', B1_idx,length(B1_scaling_arr),T1_idx,length(T1_tse_arr),T2_idx,length(T2_tse_arr),tse_idx,ETL);
        
        % collect snapshot of magnetization
        % 			Mxy_tmp = (Mx + 1i*My); len = length(Mxy_tmp);
        % 			Mxy_all = [Mxy_all Mxy_tmp(1 : len) zeros(1,5000)];
        % 			Mz_all  = [Mz_all Mz zeros(1,5000)];
        
        % --------------------------------------------------
        % pre  refocusing TE fill time
        % --------------------------------------------------
        if (TE_fill_timePre(tse_idx) > 0)
            % 		fprintf('Adding %3.0f [us] for TE fill time (pre refocusing)\n',TE_fill_timePre(tse_idx)*1e6);
            B_eff_rot = B_rot;
            st = set_evolve_M_CPP_struct(M,B_eff_rot,TE_fill_timePre(tse_idx),dB0z,0,sample_z_axis,RH_flag,Inhomo_Flag,Relax_Flag,T1_tse,T2_tse,context);
            [Mx,My,Mz] = evolve_M_CPP(st);		M = [transpose(Mx),transpose(My),transpose(Mz)];
        end;
        
        % --------------------------------------------------
        % pre  refocusing RO pre-phasing & crusher gradients
        % (skipped on first pass)
        % --------------------------------------------------
        if (tse_idx > 1)
            B_eff_rot = B_rot;
            st = set_evolve_M_CPP_struct(M,B_eff_rot,Tpr,dB0z,Gpr1 + Gcrush,sample_z_axis,RH_flag,Inhomo_Flag,Relax_Flag,T1_tse,T2_tse,context);
            [Mx,My,Mz] = evolve_M_CPP(st);
            M = [transpose(Mx),transpose(My),transpose(Mz)];
            plot_M_n_Phase(0,sample_z_axis,M_init,M,sprintf('M after %2.0f pre REF pre-phasing grad',tse_idx),sprintf('Phase after %2.0f pre REF pre-phasing grad',tse_idx),{'\phi(z)'});
        end;
        
        % --------------------------------------------------
        % Refocusing RF
        % --------------------------------------------------
        if (SE_flag)
            if ((tse_idx == ceil(ETL/2)) && HyperEcho_flag)
                % this is the mid-point iteration, apply a pure 180
                B_eff_rot = B_refoc_180 + (transpose(ones(1,N_refoc))) * B_rot;
                % for next iterations (second half of the echo-train): switch the refocusing pulses' phase and flip-angle
                tmp           = B_refoc;
                B_refoc       = B_refoc_minus;
                B_refoc_minus = tmp;
                % 			Grefoc        = -Grefoc;        % not required for HyperEcho !!!
                % 			Gcrush0       = -Gcrush0;       % not required for HyperEcho !!!
            else
                if CPMG_flag  &&  (mod(tse_idx,2) == 0)  &&  ~HyperEcho_flag
                    % flip refocusing RF on even echoes -- THIS IS IMPORTANT IN ORDER TO OBSERVE THE HIGHER-2ND-ECHO ARTIFACT
                    B_eff_rot = B_refoc_minus + (transpose(ones(1,N_refoc))) * B_rot;
                else
                    % Implementation of the CPMG condition but with exc-phase that is 90-deg off the refocusing phase
                    B_eff_rot = B_refoc       + (transpose(ones(1,N_refoc))) * B_rot;
                end;
            end;
            
            if (isVERSEref)
                Gradient = GverseRef;
            else
                Gradient = Grefoc*ones(1,N_refoc);
            end;
            
            st = set_evolve_M_n_acq1Dgrad_t_CPP_struct(M,B_eff_rot,dte_refoc,dB0z,Gradient,sample_z_axis,RH_flag,0,Inhomo_Flag,Relax_Flag,T1_tse,T2_tse,context);
            [Mx,My,Mz,dummy1,dummy2] = evolve_M_n_acq1Dgrad_t_CPP(st);
            M = [transpose(Mx),transpose(My),transpose(Mz)];
            plot_M_n_Phase(0,sample_z_axis,M_init,M,sprintf('M after %2.0f REF RF',tse_idx),sprintf('Phase after %2.0f REF RF',tse_idx),{'\phi(z)'});
        end;
        
        % --------------------------------------------------
        % post refocusing RO pre-phasing & crusher gradients
        % --------------------------------------------------
        B_eff_rot = B_rot;
        st = set_evolve_M_CPP_struct(M,B_eff_rot,Tpr,dB0z,Gpr2 + Gcrush,sample_z_axis,RH_flag,Inhomo_Flag,Relax_Flag,T1_tse,T2_tse,context);
        [Mx,My,Mz] = evolve_M_CPP(st);
        M = [transpose(Mx),transpose(My),transpose(Mz)];
        plot_M_n_Phase(0,sample_z_axis,M_init,M,sprintf('M after %2.0f post REF pre-phasing grad',tse_idx),sprintf('Phase after %2.0f post REF pre-phasing grad',tse_idx),{'\phi(z)'});
        
        % --------------------------------------------------
        % post refocusing TE fill time
        % --------------------------------------------------
        if (TE_fill_timePst(tse_idx) > 0)
            % 		fprintf('Adding %3.0f [us] for TE fill time (post refocusing)\n',TE_fill_timePst(tse_idx)*1e6);
            B_eff_rot = B_rot;
            st = set_evolve_M_CPP_struct(M,B_eff_rot,TE_fill_timePst(tse_idx),dB0z,0,sample_z_axis,RH_flag,Inhomo_Flag,Relax_Flag,T1_tse,T2_tse,context);
            [Mx,My,Mz] = evolve_M_CPP(st);		M = [transpose(Mx),transpose(My),transpose(Mz)];
        end;
        
        % --------------------------------------------------
        % Acquisition
        % --------------------------------------------------
        B_eff_rot = transpose(ones(1,length(ta))) * B_rot;
        acq_flag = 1;
        
        acq_grad = Ga * ones(1,length(ta));
        st = set_evolve_M_n_acq1Dgrad_t_CPP_struct(M,B_eff_rot,dta,dB0z,acq_grad,sample_z_axis,RH_flag,acq_flag,...
            Inhomo_Flag,Relax_Flag,T1_tse,T2_tse,context);
        [Mx,My,Mz,sig_real,sig_imag] = evolve_M_n_acq1Dgrad_t_CPP(st);
        M = [transpose(Mx),transpose(My),transpose(Mz)];
        sig_tmp = transpose((transpose(sig_real) + 1i*transpose(sig_imag)));
        % 	sig_tmp = zero_pad_1D(sig_tmp); check this if enabling code. when using low resolution imaging it distorts the image
        % 			sig_cell{tse_idx} = sig_tmp;
        % 			Sig = [Sig sig_tmp];
        im_tmp  = fftshift(fft(fftshift(sig_tmp)));  % don't use abs - keep the image phase... for now
%         im_cell{idx_ems ,tse_idx}  = im_tmp;
        % 			Im  = [Im  im_tmp(1 : length(im_tmp))];
        plot_M_n_Phase(0,sample_z_axis,M_init,M,sprintf('M after %2.0f ACQ',tse_idx),sprintf('Phase after %2.0f ACQ',tse_idx),{'\phi(z)'});
        
        % --------------------------------------------------
        % post acquisition rewinding
        % --------------------------------------------------
        if (Grew ~= 0)
            B_eff_rot = B_rot;
            st = set_evolve_M_CPP_struct(M,B_eff_rot,Tpr,dB0z,Grew,sample_z_axis,RH_flag,Inhomo_Flag,Relax_Flag,T1_tse,T2_tse,context);
            [Mx,My,Mz] = evolve_M_CPP(st);
            M = [transpose(Mx),transpose(My),transpose(Mz)];
            plot_M_n_Phase(0,sample_z_axis,M_init,M,sprintf('M after %2.0f post acquisition rewinding',tse_idx),sprintf('Phase after %2.0f post acquisition rewinding',tse_idx),{'\phi(z)'});
        end;
        

    % =========================================================================
    % Calcualte accumulated diffusion decay
    % =========================================================================
    
    
    
    % =========================================================================
    % Calcualte current EMC (image processing stage)
    % =========================================================================
    % 1. Concatinate all imaged one after another in a single vector -- Just for easy viewing
    % 2. Estimate the echo decay curve, modulated with stimulated and indirect echoes ('echo_amp' array)
%     integrated_loc = [];
%     echo_amp       = zeros(1,ETL);
%     Im            = zeros(length(TE_arr),length(im_tmp)*ETL);
%     for idx = 1:ETL
%         im_tmp = im_cell{idx_ems ,idx};
        if (EMC_Abs_flag)
            im_tmp = abs(im_tmp);
        end;
        im_tmp = interp1(1:length(im_tmp),im_tmp,linspace(1,length(im_tmp),length(im_tmp)*100));
%         	Im(idx_ems,:) = [Im(idx_ems,:) abs(im_tmp)];
        
        echo_amp(tse_idx) = sum(im_tmp(abs(im_tmp) > (sliceTH*max(abs(im_tmp))))) / length(abs(im_tmp));   % sum(im_tmp... w/o abs to keep the phase
        % 	echo_amp(idx) = max(im_tmp);                                                                   % works very bad...
        a = abs(im_tmp);
        a(a> sliceTH*max(a)) = max(a);
        a(a<=sliceTH*max(a)) = 0;
        integrated_loc = [integrated_loc a];
%     end;
        
        
    end; % TSE loop
    
    
    
    Integrated_loc(idx_ems,:) = integrated_loc;
    % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    %                                                        End  T2 & B1 Loop
    echo_train_modulation(idx_ems,:) = echo_amp;
    %             end; % for T2
    %         end; % for T1
    %     end; % for B1
    % end; % for B0
end, % Parfor idx_ems
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
echo_train_modulation = reshape(echo_train_modulation, T2_N,T1_N,B1_N,B0_N,ETL);
echo_train_modulation = permute(echo_train_modulation,[4 3 2 1 5]);


% -------------------
% Interpolate the EMC
% -------------------
% if (length(T2_tse_arr) > 1) && (EMC_interp > 1)
% 	[iecho_train_modulation, iT2_tse_arr] = interp_EMC(echo_train_modulation, B1_scaling_arr, T2_tse_arr, EMC_interp);
% end;

% --------------------------------------------------
% Save curves
% --------------------------------------------------
TSE_fn = sprintf('%s_%1.0fum_L_%1.0fmm',TSE_fn,required_res*1e4,Lhalf*2*10);
if (TSE_save_flag)
    save(TSE_fn);
end;

% % --------------------------------------------------
% % Save separate SE DICOMs - for simulating a single-echo SE pulse sequence for a single T2 value
% % --------------------------------------------------
% TSE_fn = sprintf('%s_%1.0fmm_L_%1.0fmm',TSE_fn,required_res*10,Lhalf*2*10);
% for idx = 1:length(TE_arr)
% 	TSE_fn1 = [TSE_fn sprintf('_TE%1.0fms',TE_arr(idx)*1e3)];
% 	a = im_cell{idx};
% 	a=a/2;
% 	dicomwrite(a,[TSE_fn '.scm']);
% end;

% =========================================================================
% Plot
% =========================================================================
if (length(T2_tse_arr) == 1) && (length(B1_scaling_arr) == 1) && (length(B0_scaling_arr) == 1)
% 	ax = 10*linspace(-Lhalf,+Lhalf,length(Im));
% 	figure;
% 	subplot(221); plot(abs(Mxy_all),'b-'); title(sprintf('Mxy-all refocuse angle = %3.1f; Refname=%s;   Gcrush=%3.1f',refoc_angle_tse(B1_idx),Refname,Gcrush)); grid;
% 	subplot(223); plot(Mz_all      ,'k-'); title(sprintf('Mz-all  refocuse angle = %3.1f; Refname=%s;   Gcrush=%3.1f',refoc_angle_tse(B1_idx),Refname,Gcrush)); grid;
% 	subplot(222); plot(abs(Sig)    ,'.-'); title(sprintf('Sig-all        angle = %3.1f;   Refname=%s;   Gcrush=%3.1f',refoc_angle_tse(B1_idx),Refname,Gcrush)); grid;
% 	subplot(224); plot(ax,abs(Im ) ,'.-'); title(sprintf('Im    refocuse angle = %3.1f;   Refname=%s;   Gcrush=%3.1f',refoc_angle_tse(B1_idx),Refname,Gcrush)); grid; xlabel('mm');

	echo_amp_T2a = get_fitted_T2(TE_arr ,echo_amp);
	echo_amp_T2b = get_fitted_T2(TE_arr(2:end) ,echo_amp(2:end));
	echo_amp     = echo_amp / echo_amp(1);

	figure;
	subplot(211); hold on;
				  plot(abs(Im            ),'k.-');   title(sprintf('Im  refocuse angle = %3.1f;  Refname=%s;   Gcrush=%3.1f',refoc_angle,Refname,Gcrush));
				  plot(abs(integrated_loc),'m-');    title(sprintf('integrated locations'));
	subplot(212); plot(abs(echo_amp      ),'b.-');   title(sprintf('cummulative integrals (real T2 = %3.0f [ms]  fitted T2 = %3.0f [ms], w/o 1st echo = %3.0f',T2*1e3,echo_amp_T2a*1e3,echo_amp_T2b*1e3));

	if (exist('fh1') == 1)
		if (~fh1) fh1 = figure; end;
		figure(fh1); hold on;
		plot(abs(echo_amp),styles{reqres_idx});      title(sprintf('cummulative integrals (real T2 = %3.0f [ms]  fitted T2 = %3.0f [ms], w/o 1st echo = %3.0f',T2*1e3,echo_amp_T2a*1e3,echo_amp_T2b*1e3));
% 		plot(abs(echo_amp),styles{refoc_angle_idx}); title(sprintf('cummulative integrals (real T2 = %3.0f [ms]  fitted T2 = %3.0f [ms], w/o 1st echo = %3.0f',T2*1e3,echo_amp_T2a*1e3,echo_amp_T2b*1e3));
		axis([0.5 ETL+0.5 0.0 1.3]);
		legend(lgd1);
	end;
else
	% Disable this plot for now -- it's excessively time-consuming
	figure; hold on; lgd = {};
	for idx0 = 1:1:length(B0_scaling_arr)
		for idx1 = 1:1:length(B1_scaling_arr)
			for idx2 = 1:1:length(T1_tse_arr)
				for idx3 = 1:1:length(T2_tse_arr)
					st_idx = mod((idx0+idx1+idx2+idx3)-2,length(styles)-1)+1;
					emc_vec = abs(squeeze(echo_train_modulation(idx0,idx1,idx2,idx3,:)));
					emc_vec = emc_vec / emc_vec(1);
					plot(emc_vec,styles{st_idx});
% 					plot(emc_vec);
					lgd{end+1} = ['T2=' num2str(T2_tse_arr(idx3)*1000) ' T1=' num2str(T1_tse_arr(idx2)*1000) ' B1=' num2str(B1_scaling_arr(idx1)) ' B0=' num2str(B0_scaling_arr(idx0))];
				end;
			end;
		end;
	end;
	legend(lgd);
	title(sprintf('Echo train modulation (T1 = %1.0f [ms])',T1*1000));
end;

% profile off
% profile report
time = toc;
hr   = floor(time/3600);
mn   = floor((time - hr*3600)/60);
sc   = round( time - hr*3600 - mn*60);
fprintf('\n ===== Total run time = %1.0f:%1.0f:%1.0f hr:mn:sc\n\n',hr,mn,sc);





% % % % % % ------------------------------------
% % % % % % Add noise
% % % % % % ------------------------------------
% % % % % orig_sig            = transpose(Sig);
% % % % % interp_flag         = 0;
% % % % % if (0)
% % % % % 	declare_stage('Adding noise');
% % % % % 	noise_std           = 1;             % [%]
% % % % % 	signal_lvl          = 1;
% % % % % 	complex_noise_f     = 1;
% % % % % 	Sig = add_random_noise_to_1D_sig(orig_sig,fb,noise_std,signal_lvl,complex_noise_f,reset_random_seed_f,0);
% % % % % end;
% % % % %
% % % % %
% % % % %
% % % % % % --------------
% % % % % %  plot the PSF
% % % % % % --------------
% % % % % z   = sample_z_axis;
% % % % % Lz  = 2*Lhalf;
% % % % % dza = Lz/N_acq_pts;
% % % % % PSF = (exp(1i*pi*z/Lz) .* sin(pi*z/dza)) ./ (sin(pi*z/Lz));
% % % % % acq_grid = linspace(za_axis(1),za_axis(end),N_acq_pts);
% % % % %
% % % % % figure; hold on;
% % % % % plot(z,abs(PSF),'b-');
% % % % % plot(acq_grid,zeros(1,N_acq_pts),'m^-');
% % % % % legend({'PSF','acquisition grid'});
