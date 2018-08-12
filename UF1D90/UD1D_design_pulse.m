
function [chirp_rot,phi_chirp_rot,alpha0_n,alpha1_n,alpha2] = UD1D_design_pulse(context)
set_context;

declare_stage('Design PI/2-Chirp pulse');
if (DesignPulse_Flag)
	if (Inhomo_Fix_Flag == 0)
		P_DNu0_temp = P_DNu0;
		P_DNu0      = P_DNu0*0;
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
		P_DNu0 = P_DNu0_temp;
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
		error('WRONG PARAMETERS PASSED TO SLR FUNCTION');
		chirp_rot_SLR = convert_pulse_to_SLR(z_axis,Tp,Phi_e_sq,OmegaE,ze_of_OmegaE,chirp_rot,sw,tpwr90,pw90,RF_ESP_SLR,pulse_fname);
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

		for idx = 1:Nsp
			% Calculate the pulse frequency range
			Oi_n(idx) = Oi_Hz + (idx-1)*dO_Hz;                      % [Hz]
			Of_Hz = Oi_n(idx) + dO_Hz;                              % [Hz]
			alpha0_n(idx) = - Tp * ((Oi_n(idx))^2) / (2*dO_Hz);
			alpha1_n(idx) = gammaHz * Ge(3) * Tp * Of_Hz / dO_Hz;
		end;
		alpha2 = - ((gammaHz*Ge(3))^2) / (2*rate);

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
				ax = z_axis(n(sp_idx)+1 : n(sp_idx+1));                                % Spatial range of current SP
				mid_val = mean([z_axis(n(sp_idx)+1), z_axis(n(sp_idx+1))]);            % Middle value of the range
				Phi_e_sq_rad((end+1):n(sp_idx+1)) = 2*pi*alpha2*((ax - mid_val) .^ 2); % Frequency response
			end;
		else
			Phi_e_sq_rad = 2*pi*alpha2*(z_axis.^2);
		end;

		if (DEBUG_FLAG > 2)
			figure;
			subplot(3,1,1); plot(z_axis      ,OmegaE*1E-3 ,'.-'); title('OmegaE(z)  '); xlabel('z-axis [cm]'   ); ylabel('Omega(z) [kHz]'       ); grid;
			subplot(3,1,2); plot(OmegaE*1E-3 ,z_of_OmegaE ,'.-'); title('z(OmegaE)  '); xlabel('Omega(z) [kHz]'); ylabel('z-axis [cm]'          ); grid;
			subplot(3,1,3); plot(z_axis      ,Phi_e_sq_rad,'.-'); title('Phi e sq(z)'); xlabel('z-axis [cm]'   ); ylabel('Phi_e_sq_rad(z) [rad]'); grid;
		end;
		[chirp_rot,rfpwr,rffpwr] = convert_pulse_to_SLR(z_axis,Tp,rfwdth,Phi_e_sq_rad,OmegaE,OmegaE_axis,...
		                                                z_of_OmegaE,swHz,isPI,tpwr90,pw90,RF_ESP_SLR,DEBUG_FLAG);
		phi_chirp_rot = phase(chirp_rot);
	else % if (SLR_flag)
% 		Of = Oi_Hz + rate*te;
% 		Of = Of_Hz(end);
% 		DO = Of_Hz-Oi_Hz;
		fprintf('DO=%3.3f [kHz];  R=%3.3f [kHz/ms];  Tp=%3.3f [ms]\n',DO_Hz*1e-3,rate*1e-6,(te(end)+dte)*1e+3);
		[chirp_rot, phi_chirp_rot] = gen_Chirp_pulse(te,Oi_Hz,rate,B1ampCoeff,context);
		alpha0_n = 0;
		alpha1_n = 0;
		alpha2   = 0;

		% Windowing
		if (RFnGRAD_Win_Flag)
			chirp_rot = window_rf_WURST(chirp_rot,te,RF_ESP,0);%DEBUG_FLAG);
		end;
	end;
end;

return;
