function [chirp90,rf90pwr,rf90fpwr,dO,R90,alpha0,alpha1,alpha2] = gen_chirp_90_pulse(...
                                Zie,Zfe,GPEe,Tp,dte,rf90wdth,SLR_flag,tpwr90,pw90,RF_ESP)
set_globals;

Oi  = gammaHz*GPEe*Zie;
Of  = gammaHz*GPEe*Zfe;
dO  = Of - Oi;
R90 = dO / Tp;

te = 0:dte:(rf90wdth - dte);
z_axis = linspace(Zie,Zfe,length(te));

alpha0  = - Tp * (Oi^2) / (2*dO);                   % [none]
alpha1  = gammaHz * GPEe * Tp * Of / dO;            % [1/cm]
alpha2  = - ((gammaHz*GPEe)^2) / (2*R90);           % [1/cm^2]

if (SLR_flag)
	OmegaE       = gammaHz*GPEe*z_axis;
	OmegaE_axis  = OmegaE;
	z_of_OmegaE  = interp1(OmegaE,z_axis,OmegaE);
	Phi_e_sq_rad = 2*pi*alpha2*(z_axis.^2);
	swHz         = 1/dte;
	isPI         = 0;
	if (DEBUG_FLAG >=2)
	figure;
	subplot(3,1,1); plot(z_axis     ,OmegaE*1E-3 ,'.-'); title('OmegaE(z)  '); xlabel('z-axis [cm]'   ); ylabel('Omega(z) [kHz]'); grid;
	subplot(3,1,2); plot(OmegaE*1E-3,z_of_OmegaE ,'.-'); title('z(OmegaE)  '); xlabel('Omega(z) [kHz]'); ylabel('z-axis [cm]'   ); grid;
	subplot(3,1,3); plot(z_axis     ,Phi_e_sq_rad,'.-'); title('Phi e sq(z)'); xlabel('z-axis [cm]'   ); ylabel('[rad]'         ); grid;
	end;
	disp(sprintf('Designing Chirp 90 SLR pulse. Tp=%3.3f, rf90wdth=%3.3f',Tp,rf90wdth));
	[chirp90,rf90pwr,rf90fpwr] = convert_pulse_to_SLR(z_axis,Tp,rf90wdth,Phi_e_sq_rad,...
	                                              OmegaE,OmegaE_axis,z_of_OmegaE,swHz,...
	                                              isPI,tpwr90,pw90,RF_ESP,DEBUG_FLAG);
else
	[chirp90, phi_chirp90] = gen_Chirp_pulse(te,Oi,R90,B1Chirp90Coeff);

	% Windowing
	chirp90 = window_rf(chirp90,te,2.1,RF_ESP,DEBUG_FLAG,'lala');

	% Calibration
	[rf90pwr,rf90fpwr] = calib_chirp_pulse(tpwr90,pw90,R90,0);
end;

return;