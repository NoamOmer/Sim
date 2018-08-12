
function MultSP_1DUF_pp(FIDnum)
% clear all;
clcl;
FIDnum = 592;
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
dZa = sqrt((Lz/Nsp)/(gammaHz*Ge*Tp));              % [cm]
dt = 1E-3/sw;                                      % [sec]

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
[sig_imag_filtered, phase_sign] = MultSP_reconstruct_profile(sig,Nsp,np,Ga,Ta,Ge,Tp,Lz,alpha0_n,alpha1_n,alpha2,...
                                                             ppWinExpVal,ppWinESP);

% -----------------------------------------------------------------------
% Reconstruct interface region - Op #2 - iteratively fill in each SP side
% according to the predicted phase, on the other SP side
% -----------------------------------------------------------------------
% sig_imag_filtered = sig_imag_filtered/max(sig_imag_filtered);
s0_orig = sig_imag_filtered;
s0      = sig_imag_filtered;

% Interpolate 's0' to be able to achieve higher sampling rate -- this does not
% work since the np should be set to its maximal value to start with, not during PP.
interp_flag = 0;
if (interp_flag)
	Npts = length(s0)*8;
	if (mod(Npts,2) ~= 0)
		Npts = Npts + 1;
	end;
	s0 = interp1(1:length(s0),s0,linspace(1,length(s0),Npts));
% 	s0 = [s0(7:end) zeros(1,6)];  % shift the infetface region to the center of s0
else
	Npts = 2*np;
end;

% % % Create a reference sample
% Ones
s1 = ones(1,length(s0));             % Reference signal - compensation, calculated on a 'demi' reference object

% %  Smoothed rect
% s1_z_axis = linspace(-round(length(s0)/2),+round(length(s0)/2),length(s0));
% s1_dz = 2.6/1603;
% s1_z_axis = (-2 : s1_dz : 2+s1_dz);
% s1  = exp(-70E+47 * ((s1_z_axis*s1_dz).^18));
% s1 = interp1(1:length(s1),s1,linspace(1,length(s1),Npts));

s1_orig = s1;
figure; plot(s1,'.-'); title('Reference sample');

nMz = Npts/2;                             % Number of complex points
m   = 1:2*nMz;
dta = Ta/nMz;
z_m  = -Lz/2 + m*Lz/(2*nMz);
dz   = abs(z_m(2) - z_m(1));

for iter_idx = 1:1                        % Iterations loop

	% Initialize the demi (reference) profile
	if (1)  % same profile each iteration
		demi_s = s1_orig;
	else    % incrementing profile
		demi_s = abs(s1)/max(s1);
	end;
	% if (iter_idx == 2)
	% 	demi_s(100:114) = linspace(demi_s(101),demi_s(114),15);
	% end;
	% if (iter_idx == 3)
	% 	demi_s(102:112) = linspace(demi_s(102),demi_s(112),11);
	% end;	
	if (DEBUG_FLAG >= 3)  figure; plot(demi_s,'.-'); title(sprintf('Demi signal (iter #%1.0f)',iter_idx));  end;

	% Reset the contribution vector
	this_half_contribution_vec = 0;

	for n = 1:Nsp                             % loop over SPs

		% 1.1   Calculate the excitation phase profile, relating to the current SP
		ex_phase  = alpha2 * (z_m.^2) + alpha1_n(n) * z_m + alpha0_n(n);
	% 	ex_phase2 = alpha2 * (z_m.^2) + alpha1_n(2) * z_m + alpha0_n(2);
	% 	ex_phase3 = ex_phase1 + ex_phase2;
	% 	figure;plot(z_m,ex_phase1,'b.-',z_m,ex_phase2,'r.-',z_m,ex_phase3,'m.-'); legend({'1','2','1+2'});

		% 1.2   Calculate the sample indices that relate to current SP
		this_half_idxs = (1:nMz) + (n-1)*nMz;

		% 2     Set the predicted phase on the reconstructed signal
		%       For each acquired point, add the contribution from the "other half" of the sample

		for idx = 1:nMz                       % Loop over all point in this half of the sample

			% 2.1   Add the current acquisition phase to the initial excitation phase
			%       convert to [rad] & wrap into range [-pi .. +pi]
			sp_phase = ex_phase + gammaHz*Ga*dta*(idx-1).*z_m;		%DEBUG	figure(1001); hold on; plot(ex_phase,'b.-'); plot(sp_phase,'r.-'); legend({'ex_phase','sp_phase'});
			sp_phase = sp_phase*2*pi;
			sp_phase_wrapped = mod(sp_phase,2*pi) - pi;

			% 2.2   Set the phase profile on the reference profile
			demi_s = abs(demi_s) .* exp(phase_sign*i*sp_phase_wrapped);

			% 2.3   "Acquire" (i.e., sum up) the contribution from "this half" of the sample
			this_half_contribution = sum(demi_s(this_half_idxs))*dz;

			% 2.4   Add "this half" contribution of the current SP to the reference signal.
			%       Start by calculating the current location (i.e., index) of the SP
			%       on this side of the sample
			sp_idx = 1 + n*nMz - idx;  % scalar

			% Note: The 'this_half_contribution_vec' vector is not overwritten between iterations.
			%       Each iteration fills only its own SP indices
			this_half_contribution_vec(sp_idx) = this_half_contribution;
		end;
	end; % SP loop

	% 3   Prepare the reference contribution vector
	demi_s = abs(this_half_contribution_vec);
	% demi_s = demi_s ./ max(demi_s);

	% 6   Reconstruct the initial signal according to the reference contribution vector
	% 6.1 Correct only around the interface region - use two methods for this localization
	%     1  Window the reference vector around the interface region --> Filter HF components
	%     2  FT-->Window--> IFT the reference vector to remove HF components --> Window
	a1  = demi_s;
	la1 = length(a1);
	quart_loc = round(la1/4);
	if (la1 >= 500)
		shift = round(la1/10);
		a2 = a1 / mean(a1(quart_loc-shift:quart_loc+shift));
	else
% 		a2 = a1 / mean(a1(quart_loc));
		a2 = a1 / max(a1);
	end;
	a3 = a2 - 1;

	% Window the ref vector around the center
	b1 = window_1Dvec(a3,15+(iter_idx-1)*0,10,0,'none');
	b1 = b1 + 1;

	% Filter out the very HF components
	fft_b1 = fftshift(fft(fftshift(b1)));
	fil_fft_b1 = window_1Dvec(fft_b1,3,10,0,'none');
	b1 = ifftshift(ifft(ifftshift(fil_fft_b1)));

	if (DEBUG_FLAG >= 1)  
	figure; plot(1:la1,a1,'b.-',1:la1,a2,'r.-',1:la1,a3,'k.-',1:la1,b1,'k-',1:la1);
	legend({'Init ref vec','a2','a3','Final ref vec'}); grid; title('Filter #2.1');
	end;
	
% 	% 6.2 Scale the initial signal according to the filtered reference vector
% 	%     Multiply the result by the profile, used for the reference scan
% 	tmp1 = this_half_contribution_vec(1:107) + this_half_contribution_vec(108:214);
% 	[tmp2, phase_sign] = MultSP_reconstruct_profile(tmp1,Nsp,np,Ga,Ta,Ge,Tp,Lz,alpha0_n,alpha1_n,alpha2,...
% 													ppWinExpVal,ppWinESP);
% 	figure; plot(abs(tmp2),'.-'); title('Ref vec --> merged --> separated');
% 
% 	b2 = tmp2 / max(tmp2);
% 	b2 = b2 - 1;
% 	b2 = window_1Dvec(b2,60,10,1,'none');
% 	b2 = b2 + 1;
% 	s2 = abs(s0_orig) / max(abs(s0_orig));
% 	s2 = s2 ./ abs(b2);
% 	figure;
% 	plot(s2/max(s2)                    ,'r.-'); hold on;
% 	plot(abs(s0_orig)/max(abs(s0_orig)),'b.-'); legend({'Reconstructed s2','Original S'});
% 	title('Reconstructed S2');

	s1   = (abs(s0)./(abs(b1).^1));

	% % 6.3 Remove redundant high frequency components
	% fft_s1 = fftshift(fft(fftshift(s1)));
	% fil_fft_s1 = window_1Dvec(fft_s1,20,5,0,'none');
	% s1 = abs(ifftshift(ifft(ifftshift(fil_fft_s1))));
	% s1 = s1 / max(s1);

	% figure;
	% subplot(2,1,2); plot(s1/max(s1)          ,'r.-'); legend({'Reconstructed S (Op #2.1)'}); hold on;
	% subplot(2,1,1); plot(abs(s0)/max(abs(s0)),'b.-'); legend({'Original S'});
	% title('Reconstructed S (Op #2.1)');
	if (DEBUG_FLAG >= 2)
	figure;
	plot(s1/max(s1)                    ,'r-'); hold on;
	plot(s0/max(s0)                    ,'m-');
	plot(abs(s0_orig)/max(abs(s0_orig)),'b-'); %legend({'Reconstructed S','Previous iteration S','Original S'});
	title(sprintf('Reconstructed S (Op #2.1) (iter %1.0d)',iter_idx)); axis([1 length(s0) 0 1]); grid;
	end;
	
	% Smooth the interface region. Smoothed region width is ~ 1 pixel
	s_tmp = s1/max(s1);
	mid_point = round(length(s_tmp)/2);
	pixel_len = Lz / length(s_tmp);
	shift     = round(0.5*dZa/pixel_len);
	s_tmp((mid_point-shift) : (mid_point+shift)) = smooth(s_tmp((mid_point-shift) : (mid_point+shift)),shift);
	s_tmp = s_tmp / max(s_tmp);
	figure;
	plot(s_tmp                         ,'r-'); hold on;
	plot(s0/max(s0)                    ,'m-');
	plot(abs(s0_orig)/max(abs(s0_orig)),'b-'); %legend({'Reconstructed & smoothed S','Previous iteration S','Original S'});
	title(sprintf('Reconstructed & smoothed S (Op #2.1) (iter %1.0d)',iter_idx)); axis([1 length(s0) 0 1]); grid;
	s0 = s_tmp;%s1;
	% figure; plot(s1/max(s1)          ,'r.-'); title('Reconstructed S (Op #2.1)');
	% figure; plot(abs(s0)/max(abs(s0)),'b.-'); title('Original S');

end; % # of iterations loop

return;


% s2 = abs(s0)./abs(b2);

% % 	else % Why does this method give better result??? it shouldn't - the 1st method should be better...
% 	a1  = abs(this_half_contribution_vec);
% 	la1 = length(a1);
% 	quart_loc = round(la1/4);
% 	a2 = a1 / mean(a1(quart_loc-100:quart_loc+100));
% 	a3 = a2 - 1;
% 	b2 = window_1Dvec(a3,7,15,0,'none');
% 	b2 = b2 + 1;
% 	fft_b2 = fftshift(fft(fftshift(b2)));
% 	fil_fft_b2 = window_1Dvec(fft_b2,3,10,0,'none');
% 	b2 = ifftshift(ifft(ifftshift(fil_fft_b2)));
% 	figure; plot(1:la1,a1,'b-',1:la1,b2,'k-'); legend({'Init ref vec','Filt. ref. vec'}); grid;
% 	title('Filter #2.2');
% 	end;


% figure;
% subplot(2,1,2); plot(s2/max(s2)          ,'r.-'); legend({'Reconstructed S (Op #2.2)'}); hold on;
% subplot(2,1,1); plot(abs(s0)/max(abs(s0)),'b.-'); legend({'Original S'});
% % legend({'Reconstructed S (Op #2.2)','Original S'});
% title('Reconstructed S (Op #2.2)');
