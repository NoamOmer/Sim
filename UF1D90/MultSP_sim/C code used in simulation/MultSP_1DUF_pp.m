
function MultSP_1DUF_pp(FIDnum)
clear all;
% clcl;
FIDnum = 587;
set_globals;
ppWinExpVal = 15;
ppWinESP    = 3;

fn = sprintf('%d.fid',FIDnum);
if (exist(fn,'file') ~= 2)
	errordlg(sprintf('File %s does not exist',fn),'File Error');
	return;
end;

run(['set_' fn(1:end-4)]);

% -----------------------------------------------
% Initialize parameters
% -----------------------------------------------
dt      = 1E-3/sw;                                 % [sec]

Oi = Ge*Zie*gammaHz;                               % [Hz]
dO = rf_R*Tp;                                      % [Hz]
for idx = 1:Nsp
	% Calculate the pulse frequency range
	Oi_n(idx) = Oi + (idx-1)*dO;                   % [Hz]
	Of = Oi_n(idx) + dO;                           % [Hz]
    
    alpha0_n(idx) = - Tp * ((Oi_n(idx))^2) / (2*dO);
	alpha1_n(idx) = gammaHz * Ge * Tp * Of / dO;
end;
alpha2 = - ((gammaHz*Ge)^2) / (2*rf_R);

fid = load(fn);
sig = transpose(fid(:,1) + i*(fid(:,2)));

% Reverse sig to match the excitation direction
sig = flipdim(sig,2);

if (DEBUG_FLAG >= 3)
figure; plot(abs(sig)  ,'.-'); title(sprintf('Raw %d.FID (Amp)'  ,FIDnum)); xlabel('z-axis'); set_gca;
figure; plot(phase(sig),'.-'); title(sprintf('Raw %d.FID (phase)',FIDnum)); xlabel('z-axis'); set_gca;
end;

% -----------------------------------------------
% Post process
% -----------------------------------------------
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

	% 2.3   Filter: leave only the DC part
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
% 	figure;
% 	plot(k, abs(fft_sig_filtered)   ,'.-'); 
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

% ----------------------------------------------------------------------
% Reconstruct interface region - Op #1 - according to excitation profile
% ----------------------------------------------------------------------
if (1)
	s0 = sig_imag_filtered;
	nMz = np;                               % Number of complex points
	m   = 1:2*nMz;
	z_m  = -Lz/2 + m*Lz/(2*nMz);

	fn = sprintf('exc_amp_comp_%3.0f.mat',FIDnum);
	load(fn);

	sample_Zi        = +2.00*sign(Zie);     % [cm]      Sample z-axis
	sample_Zf        = +2.00*sign(Zfe);
	sample_Lz        = abs(sample_Zi-sample_Zf);
	margin = round( (1 - Lz / sample_Lz) * length(xd) / 2 );
	sig_xd = xd((margin+1) : (end-margin));
	sig_yd = yd((margin+1) : (end-margin));
	xd_interp = linspace(sig_xd(1),sig_xd(end),length(sig_filtered));
	yd_interp = interp1(sig_xd,sig_yd,xd_interp);

 	s_op1 = s0 ./ yd_interp;
	figure; hold on;
	plot(z_m,abs(s0),'b.-',z_m,abs(s_op1),'r.-');
	title('Interpolation result (Op #1)'); legend({'s_0','s final'});
end;

% -----------------------------------------------------------------------
% Reconstruct interface region - Op #2 - iteratively fill in each SP side
% according to the predicted phase, on the other SP side
% -----------------------------------------------------------------------
if (1)
% sig_imag_filtered = sig_imag_filtered/max(sig_imag_filtered);
s0 = sig_imag_filtered;

% Interpolate s to be able to achieve higher sampling rate -- this does not
% work since the np should be set to its maximal value to start with and
% not during post-processing.
% interp_flag = 0;
% if (interp_flag)
% 	Npts = length(sig_xd);
% 	if (mod(Npts,2) ~= 0)
% 		Npts = Npts + 1;
% 	end;
% 	s0 = interp1(1:length(s0),s0,linspace(1,length(s0),Npts));
% % 	s0 = [s0(7:end) zeros(1,6)];
% else
	Npts = 2*np;
% end;
incre_s = s0;                             % Incremented signal - incremented between iterations
final_s = s0;                             % Evaluated once according to the final value of 'incre_s'
demi_s  = ones(1,length(s0));             % Reference signal - compensation, calculated on a 'demi' reference object
dz_factor = 35;                           % Integration result factor
for iter_idx = 1:1                        % # of Iteration loop
	nMz = Npts/2;                         % Number of complex points
	m   = 1:2*nMz;
	dta = Ta/nMz;
	z_m  = -Lz/2 + m*Lz/(2*nMz);
	dz   = abs(z_m(2) - z_m(1));
	dz   = dz*dz_factor;

	for n = 1:Nsp                         % # of SP loop
		
		% 1   Calculate the excitation phase part of the current SP
		ex_phase = alpha2 * (z_m.^2) + alpha1_n(n) * z_m + alpha0_n(n);

		% 2   Set the predicted phase on the reconstructed signal
		%     For each acquired point, add the contribution from the "other half" of the sample
		other_half_idxs = (1:nMz) + (2-n)*nMz;
		 this_half_idxs = (1:nMz) + (n-1)*nMz;

		for idx = 1:nMz                   % Loop over all point in this half of the sample
			
			% 2.1   Add the current acquisition phase to the basic excitation phase
			sp_phase = ex_phase + gammaHz*Ga*dta*(idx-1).*z_m;
%DEBUG		figure(1001); hold on; plot(ex_phase,'b.-'); plot(sp_phase,'r.-'); legend({'ex_phase','sp_phase'});
			sp_phase = sp_phase*2*pi;
			sp_phase_wrapped = mod(sp_phase,2*pi) - pi;    % convert to [rad] & wrap into range [-pi .. +pi]

			% 2.2   Set the full phase on the reconstructed signal
			cur_s  = abs(incre_s) .* exp(phase_sign*i*sp_phase_wrapped);
			demi_s = abs(demi_s ) .* exp(phase_sign*i*sp_phase_wrapped);
			s0     = abs(s0     ) .* exp(phase_sign*i*sp_phase_wrapped);

			% 2.3   Sum up the contribution from the "other half" of the sample
%DEBUG		other_half_cum_contribution = cumsum(incre_s(other_half_idxs))*dz;
%DEBUG		figure(1002); plot(z_m(other_half_idxs),abs(other_half_cum_contribution),'.-');
			other_half_contribution = sum(cur_s(other_half_idxs))*dz;
			 this_half_contribution = sum(demi_s(this_half_idxs))*dz;

			% 2.4   Add the "other half" contribution of the current SP to the initial signal
			sp_idx = 1 + n*nMz - idx;
			
			% Note: The 'this_half_contribution_vec' vector is not overwritten between iterations.
			%       Each iteration fills only its own SP indices
			incre_s(sp_idx) = cur_s(sp_idx) + other_half_contribution;
			final_s(sp_idx) =    s0(sp_idx) + other_half_contribution;
			this_half_contribution_vec(sp_idx) = this_half_contribution;
		end;
	end;

	% 3   Set the reference contribution set
	demi_s = abs(this_half_contribution_vec);
% 	demi_s = demi_s ./ max(demi_s);
	
	% 4   Filter the incremented and final signals - just for increasing SNR
	ppWinExpVal = ppWinExpVal * 1;
	ppWinESP    = ppWinESP    * 1;
	
	fft_final_s            = fftshift(fft(fftshift(abs(final_s))));
	fft_final_s_filtered   = window_1Dvec(fft_final_s,ppWinExpVal,ppWinESP,0,'none');
	final_s_post_filtering = ifftshift(ifft(ifftshift(fft_final_s_filtered)));

	fft_incre_s            = fftshift(fft(fftshift(abs(incre_s))));
	fft_incre_s_filtered   = window_1Dvec(fft_incre_s,ppWinExpVal,ppWinESP,0,'none');
	incre_s_post_filtering = ifftshift(ifft(ifftshift(fft_incre_s_filtered)));

	% 5   Plot the incremented and final signals - w/ & w/o  filtering
	z = linspace(Zie,Zfe,length(s0));
	figure;
	subplot(2,1,1); hold on;
	plot(z,abs(s0),'b-', z,abs(final_s),'r-');%, z,abs(incre_s),'k-');
	title(sprintf('Reconstructed S (Op #2.0) (iteration #%2.0f)',iter_idx)); legend({'S_0','S final'});%,'S Incre'});
	subplot(2,1,2);
	plot(z, abs(final_s) - abs(s0),'.-');
	title('S_{final} - S_0 difference');
	
	figure;
	subplot(2,1,1); hold on;
	plot(z,abs(s0),'b-', z,abs(final_s_post_filtering),'r-', z,abs(incre_s_post_filtering),'k-');
	title(sprintf('Post Filtering Reconstructed S (Op #2.0) (iteration #%2.0f)',iter_idx)); legend({'S_0','S final','S Incre'});
	subplot(2,1,2);
	plot(z, abs(final_s_post_filtering) - abs(s0),'.-');
	title('Post Filtering S_{final} - S_0 difference');

	% 6   Reconstruct the initial signal according to the reference contribution vector
	% 6.1 Correct only around the interface region - use two methods for this localization
	%     1  Window the reference vector around the interface region
	%     2  FT-->Window--> IFT the reference vectorin order to remove high-frequency components
% 	if (0)
		a1 = abs(this_half_contribution_vec);
		la = length(a1);
		quart_loc = round(la/4);
		a2 = a1 / mean(a1(quart_loc-100:quart_loc+100));
		a3 = a2 - 1;
		b1 = window_1Dvec(a3,16,5,0,'none');
		b1 = b1 + 1;
		fh = figure; plot(1:la,a1,'b.',1:la,a2,'r.',1:la,a3,'k.',1:la,b1,'k-');
		legend({'Init ref vec','a2','a3','Final ref vec'}); grid
% 	else % Why does this method give better result??? it shouldn't - the 1st method should be better...
		a = abs(this_half_contribution_vec);
		la = length(a);
		fft_a = fftshift(fft(fftshift(a)));
		fil_fft_a = window_1Dvec(fft_a,16,5,0,'none');
		b2 = ifftshift(ifft(ifftshift(fil_fft_a)));
		figure; plot(1:la,a,'b-',1:la,b2,'k-'); legend({'Init ref vec','Filt. ref. vec'}); grid;
% 	end;

	% 6.2 Scale the initial signal accorsing to the filtered reference vector
	s1 = abs(s0)./abs(b1);
	s2 = abs(s0)./abs(b2);
	
	figure;
	plot(s1/max(s1)          ,'r.-'); hold on;
	plot(abs(s0)/max(abs(s0)),'b.-');
	legend({'Reconstructed S (Op #2.1)','Original S'});	title('Reconstructed S (Op #2.1)');

	figure;
	plot(s2/max(s2)          ,'r.-'); hold on;
	plot(abs(s0)/max(abs(s0)),'b.-');
	legend({'Reconstructed S (Op #2.2)','Original S'}); title('Reconstructed S (Op #2.2)');
end;
end; % if interpolate according to 2nd option

% -------------------------------------------------------------------------
% Interpolate interface region - Op #3 - add to each point in the interface
% region the contribution from the other sample side. Confine the process
% within a region, equal to a finite number of pixel sizes
% -------------------------------------------------------------------------
if (1)
nMz = np;                                % Number of complex points
m   = 1:2*nMz;
z_m  = -Lz/2 + m*Lz/(2*nMz);
dz   = abs(z_m(2) - z_m(1));
dZa = sqrt((Lz/Nsp)/(gammaHz*Ge*Tp));
n_pixels = 2;                              % Width of reconstructed region
n_iter = ceil((n_pixels*dZa/dz)/2);        % Extra /2 factor since each side misses only half the contribution
s0    = sig_imag_filtered;
s_op3 = sig_imag_filtered;
left_idx  = round(length(s0)/2);
right_idx = left_idx + 1;
% Weight the contribution so that regions, far from the interface will have less effect 
weight_cont = linspace(0.10,0,n_iter);             % Linear weighting
% weight_cont = exp(weight_cont.^0.1);             % Exponential weighting
% weight_cont = weight_cont ./ max(weight_cont);   % Normalize the weighting factor

for idx = 1:n_iter
	s_op3(left_idx -idx+1) = s_op3(left_idx -idx+1) + sum(s0(right_idx:+1:(right_idx+n_iter-idx)) .* weight_cont(idx:end));
	s_op3(right_idx+idx-1) = s_op3(right_idx+idx-1) + sum(s0(left_idx :-1:(left_idx -n_iter+idx)) .* weight_cont(idx:end));
end;

figure; hold on;
plot(z_m,abs(s0),'b.-',z_m,abs(s_op3),'r.-');
title('Interpolation result (Op #3)'); legend({'s_0','s final'});
end;

return;

