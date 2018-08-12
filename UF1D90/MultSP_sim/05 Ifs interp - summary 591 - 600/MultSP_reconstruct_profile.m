
function [sig_imag_filtered, phase_sign] = MultSP_reconstruct_profile(sig,Nsp,np,Ga,Ta,Ge,Tp,Lz,...
	                                                                  alpha0_n,alpha1_n,alpha2,...
																	  ppWinExpVal,ppWinESP)
set_globals;

nMz = np;                                          % Number of complex points
m   = 1:nMz;
dta = Ta/nMz;

fft_sig = fftshift(fft(fftshift(sig)));

% Loop over stationary points - start with the 1st one, corresponding to y_axis(1 -> ...)
% uiwait(msgbox('Have you set the correct value for SIMULATION parameter?'));
if (exist('SIMULATION','var') && (SIMULATION == 1))
	phase_sign = -1;
else
	phase_sign = +1;
end;

for n = 1:Nsp
	z_m = -Lz/2 + ((m + (n-1)*nMz) - 1)*Lz/(2*nMz);

	% 2.1   Perform SP specific phase cancellation on the column vector
	sp_phase = alpha2 * (z_m.^2) + alpha1_n(n) * z_m + alpha0_n(n) + gammaHz*Ga*dta*(nMz-m).*z_m;
	sp_phase = sp_phase*2*pi;
	sp_phase_wrapped = mod(sp_phase,2*pi) - pi;    % convert to [rad] & wrap into range [-pi .. +pi]
% 	sp_phase_wrapped = sp_phase;
	sig_no_sp_specific_phase = sig .* exp(phase_sign*i*sp_phase_wrapped);

	% 2.2   Fourier Transform along the PE-axis
	fft_sig_no_sp_phase = fftshift(fft(fftshift(sig_no_sp_specific_phase)));
	Dx = 0.5*Lz;
	dx = Lz / (2*nMz);
	x  = linspace(-Dx/2,Dx/2,nMz);
	Dk = 1/dx;
	dk = 1/Dx;
	k  = linspace(-Dk/2,Dk/2,nMz);
	k0_tag = gammaHz*Ge*Tp;
	k0_tag = mod(k0_tag,Dk);
	disp(sprintf('\n nSP=%1.0f:\n  Dk = [%3.1f ... %3.1f]\n  k0 = %3.1f\n',n,-Dk/2,+Dk/2,k0_tag));

	% 2.3   Filter: leave only the current SP echo (should be located around the DC)
	fft_sig_filtered = window_1Dvec(fft_sig_no_sp_phase,ppWinExpVal,ppWinESP,0,'none');
% 	fft_sig_filtered = abs(fft_sig_filtered) .* transpose(exp(i*smooth(phase(fft_sig_filtered))));

	% 2.4   Inverse FT the filtered vector
	sig_post_filtering = ifftshift(ifft(ifftshift(fft_sig_filtered)));

    if (DEBUG_FLAG >= 3)
		f1 = figure;
		subplot(2,1,1);
		plot(z_m,phase(sig),'.-',z_m,phase(sig_no_sp_specific_phase),'.-',z_m,sp_phase,'.-');
		title('Column phase: Pre & Post SP specific cancallation');
		legend({'Pre-cancallation','Post-cancallation','SP phase'},'Location','Best');
		subplot(2,1,2);
		plot(z_m,angle(sig),'.-',z_m,angle(sig_no_sp_specific_phase),'.-');
		title('Column angle: Pre & Post SP specific cancallation');
		legend({'Pre-cancallation','Post-cancallation'},'Location','Best');

		f2 = figure;
		subplot(3,1,3); plot( abs(sig),'b.-');  hold on; plot( abs(sig_no_sp_specific_phase),'k.-');
		legend({'Pre (abs) ','Post (abs) '});
		subplot(3,1,2); plot(imag(sig),'b.-');  hold on; plot(imag(sig_no_sp_specific_phase),'k.-');
		legend({'Pre (imag)','Post (imag)'});
		subplot(3,1,1); plot(real(sig),'b.-');  hold on; plot(real(sig_no_sp_specific_phase),'k.-');
		legend({'Pre (real)','Post (real)'});
		title('Column vector - Pre & Post SP specific filtering');

		f3 = figure;
		subplot(3,2,1); plot(k,  abs(fft_sig)            ,'.-'); title('Initial column spectrum - Amp.'     );
		subplot(3,2,3); plot(k,phase(fft_sig)            ,'.-'); title('Initial column spectrum - Phase'    );
		subplot(3,2,5); plot(k,angle(fft_sig)            ,'.-'); title('Initial column spectrum - Angle'    );
		subplot(3,2,2); plot(k,  abs(fft_sig_no_sp_phase),'.-'); title('Column spectrum no SP phase - Amp.' );
		subplot(3,2,4); plot(k,phase(fft_sig_no_sp_phase),'.-'); title('Column spectrum no SP phase - Phase');
		subplot(3,2,6); plot(k,angle(fft_sig_no_sp_phase),'.-'); title('Column spectrum no SP phase - Angle'); xlabel('x-axis [cm^{-1}]');

		f4 = figure;
		subplot(3,2,1); plot(k,  abs(fft_sig_no_sp_phase),'.-'); title('Column spectrum before filtering - Amp.' );
		subplot(3,2,3); plot(k,phase(fft_sig_no_sp_phase),'.-'); title('Column spectrum before filtering - Phase');
		subplot(3,2,5); plot(k,angle(fft_sig_no_sp_phase),'.-'); title('Column spectrum before filtering - Angle');
		subplot(3,2,2); plot(k,  abs(fft_sig_filtered   ),'.-'); title('Column spectrum after  filtering - Amp.' );
		subplot(3,2,4); plot(k,phase(fft_sig_filtered   ),'.-'); title('Column spectrum after  filtering - Phase');
		subplot(3,2,6); plot(k,angle(fft_sig_filtered   ),'.-'); title('Column spectrum after  filtering - Angle'); xlabel('k-axis [cm^{-1}]');

		f5 = figure;
		subplot(2,2,1); plot(  abs(sig)               ,'.-'); title('Initial column vector - Amp.'         );
		subplot(2,2,3); plot(phase(sig)               ,'.-'); title('Initial column vector - Phase'        );
		subplot(2,2,2); plot(  abs(sig_post_filtering),'.-'); title('Column vector after filtering - Amp. ');
		subplot(2,2,4); plot(phase(sig_post_filtering),'.-'); title('Column vector after filtering - Phase');
    end;
    
	filter_window = window_1Dvec(ones(1,length(fft_sig_no_sp_phase)),ppWinExpVal,ppWinESP,0,'none');
	filter_window = filter_window * max(abs(fft_sig_no_sp_phase));
 	if (DEBUG_FLAG >= 3)
	figure;
	plot(k, abs(fft_sig_no_sp_phase),'.-'); hold on;
	plot(k, filter_window           ,'r-');
 	end;
	
    % 2.5   Rebuild the signal 2D matrix -- each column is built bottom up (1st row --> last row)
	sig_filtered(m + (n-1)*nMz) = sig_post_filtering;
end;

% -----------------------------------------------
% Finalize signal by taking the absolute value
% -----------------------------------------------
sig_imag          = abs(sig);
sig_imag_filtered = abs(sig_filtered);

% -----------------------------------------------
% Plot
% -----------------------------------------------
z_axis = linspace(-Lz/2,+Lz/2,np*Nsp);
figure; plot(       sig_imag         ,'.-'); title('Pre PP image (No interface reconstruction)');  xlabel('z-axis [cm]');
figure; plot(z_axis,sig_imag_filtered,'.-'); title('Post PP image (No interface reconstruction)'); xlabel('z-axis [cm]');

return;
