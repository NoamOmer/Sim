% Ver 2 : increase_2D_Hybrid_image_resolution2
%       : Remove all FID matrix manipulations (should be done before invoking the function)
% Ver 3 : increase_2D_Hybrid_image_resolution3
%       : (1) adjust axis and directionality conventions
%       : (2) assume that (Yia == Yfe) and that (Yfa == Yie)
%       : (3) assume that the matrix sig_mat is ordered so that its lines correspond to Yie --> Yfe
%       : (4) add windowing exponential and ESP values
% Ver 4 : increase_2D_Hybrid_image_resolution4
%       : (1) Re-organized phase map axes
%       : (2) Revised debugging plots (added 'ydir','normal' etc.)
%       : increase_2D_Hybrid_image_resolution5
%       : (1) Modified inhomogeneity attenuation model -- use exp(-gamma*dB0*t) instead of phasors integral
function [sig_mat_both_SR,cond_num]=increase_2D_Hybrid_image_resolution5(sig_mat_both,...
                          Yie,Yfe,Ly,initial_dy,required_dy,alpha0,alpha1,alpha2,GePE,...
                          Tp,NSE,TaPE,Ta,GaPE,SE_flag,SRExp,SRESP,SR_ones_win_flag,...
                          SRnIter,SR_zero_init_guess,SR_fix_inhomo_flag,SR_cs,ax_spect,...
                          T2PE_factor,pix_shift_ax)
set_globals;

% ----------------
%  Set parameters
% ----------------
% Adjust the required-resolution value so that the PE axis will contain a round number of pixels
new_dy = Ly / round(Ly/required_dy);

% Set the spatial axes
% 1. initial axis, according to initial resolution
% 2. new axis, according to required resolution
% 3. HI-Res axis, used for assigning weights to each pixel in the new axis, due to internal dephasing
init_y_axis  = linspace(Yie,Yfe,round(Ly/initial_dy));
y_axis       = linspace(Yie,Yfe,round(Ly/new_dy    ));

HRF  = 100;                                                   % High-resolution factor
HR_y_axis = linspace(Yie,Yfe,length(y_axis)*HRF);             % High-resolution x-axis: HRF pts per pixel

% disp(sprintf('New pixel (%3.2f[mm]) is smaller than initial pixel (%3.2f[mm]) by a factor of %3.1f',new_dy*10,initial_dy*10,initial_dy/new_dy));

% ---------------------------------------
%  Calculate the transformation matrix A
% ---------------------------------------
N = length(y_axis);                                           %          Number of pixels on the PE axis
M = NSE;                                                      %          Number acquisition points on the PE axis

% CS parameters
% Wcs       = -1e3;
% CS_factor = 1;
te        = linspace(0,Tp,length(y_axis));
te_HR     = linspace(0,Tp,length(y_axis)*HRF);

% Excitation phase
a0 = 2*pi*alpha0;                                             % [rad]
a1 = 2*pi*alpha1;                                             % [rad/cm]
a2 = 2*pi*alpha2;                                             % [rad/cm^2]
phi_e    = (a2*(   y_axis.^2) + a1*   y_axis + a0);           % [rad]    N      Excitation phase axis
phi_e_HR = (a2*(HR_y_axis.^2) + a1*HR_y_axis + a0);           % [rad]    NxHRF  High resolution ...

% phi_e_cs    = 2*pi*Wcs*te;
% phi_e_cs_HR = 2*pi*Wcs*te_HR;

% if (SE_flag), GaPE  = -GaPE;                      end;      % (only for some multi-SP fid's: 43,...)
if (SE_flag)
    phi_e    = -phi_e;
    phi_e_HR = -phi_e_HR;
end;
% if (SE_flag), phi_e_cs = -phi_e_cs; phi_e_cs_HR = -phi_e_cs_HR; end;

taPE = TaPE*linspace(0,M-1,M);                                % [sec]    Acquisition temporal axis
k    = 2*pi * gammaHz*GaPE*taPE;                              % [rad/cm]
% k_cs = 2*pi*Wcs*taPE; % [rad]
k_HR = linspace(k(1),k(end),length(k)*HRF);

if (DEBUG_FLAG >= 3)  
figure;
subplot(3,1,1); plot(phi_e,'.-');                  title('\Phi_e');
subplot(3,1,2); plot(k,'.-');                      title('k');
subplot(3,1,3); plot(phi_e + k(20) *y_axis,'k.-'); title('\Phi_e + k(t)y'); hold on;
subplot(3,1,3); plot(phi_e + k(end)*y_axis,'r.-'); legend({'\Phi_e + k(t_{20})y','\Phi_e + k(t_{end})y'});
figure;
subplot(1,2,1); imagesc(ones(M,1)*phi_e);     title('\phi_e(y)');   set(gca,'YDir','normal');
subplot(1,2,2); imagesc(transpose(k)*y_axis); title('\phi_a(y,t)'); set(gca,'YDir','normal');
end;
% return;
A    = exp(i*(ones(M,1)*phi_e    + transpose(k)*y_axis));     % [rad]    MxN (M lines, N columns)
A_HR = exp(i*(ones(M,1)*phi_e_HR + transpose(k)*HR_y_axis));  % [rad]    Mx(N*HRF) (M lines, N*HRF columns)

% A_cs    = exp(i*(ones(M,1)*phi_e_cs    + transpose(k_cs)*ones(1,N)));
% A_cs_HR = exp(i*(ones(M,1)*phi_e_cs_HR + transpose(k_cs)*ones(1,N*HRF)));

% A    = A .* CS_factor .* A_cs;
% A_HR = A_HR .* CS_factor .* A_cs_HR;


% ----------------------------------------
%  Perform the high-resolution evaluation
% ----------------------------------------
% (1)  Start by calculating the weights around the stationary point, at each time-point
%      This is done only once, and applied to all PE lines
W_f = calculate_spatio_temporal_weights(A_HR,N,M,HRF,Ta,T2PE_factor,0,SRExp,SRESP,SR_ones_win_flag);

if (DEBUG_FLAG >= 3)
% ta_tmp = linspace(0,Ta,M);

A_tmp = ones(M,1)*phi_e_HR + transpose(k)*HR_y_axis;
figure;
subplot(1,2,1); imagesc(A_tmp); set(gca,'YDir','normal'); title('A [rad]');
subplot(1,2,2); imagesc(W_f);   set(gca,'YDir','normal'); title('W f [normalized]');

figure; imagesc(W_f); set(gca,'YDir','normal');
end;
% A_LR = ones(M,1)*phi_e + transpose(k)*y_axis;
% figure;
% imagesc(A_LR.*W_f); set(gca,'YDir','normal'); title('A (weighted)');

% figure;
% subplot(1,2,1); imagesc(ones(M,1)*phi_e_HR    + transpose(k)*HR_y_axis + ...
%                         ones(M,1)*phi_e_cs_HR + transpose(k_cs)*ones(1,N*HRF)); title('A [rad] w/ CS');


if (~SR_fix_inhomo_flag)
	A  = A .* W_f;
% 	OPTS.disp = 0;
% 	eig_vals = eigs(A,2,'BE',OPTS);
% 	eig_vals1 = abs(eig(A));
% 	figure;plot(eig_vals1,'.-');
else
	A0 = A .* W_f;
end;
% max_eig  = eig_vals(1);
% min_eig  = eig_vals(2);
% cond_num = abs(max_eig / min_eig);
% cond_num = max(eig_vals1) / min(eig_vals1);
% cond_num = cond(A); % which is:  norm(A) * norm(pinv(A))
cond_num = 0;

% (2) Loop over all RO points and apply the SR algorithm
check_idx = round(size(sig_mat_both,1)/2); % 75;
for RO_idx = 1:size(sig_mat_both,1)
	plot_f = 0;
	if (RO_idx == check_idx)
		disp(''); plot_f = 1;
	end;

	% (2.1)  Extract an initial guess
	% - This might be done by FT-ing S, filtering all spectral components that are above
	%   the current resolution (meanining higher than 1/initial_dx) and IFT-ing the result.
	S  = transpose(sig_mat_both(RO_idx,:));                             % M
    
    % Try to avoid ringing artifacts from the sample edge in case our FOV
    % is not larger than the sample.
%     dbstop in increase_2D_Hybrid_image_resolution5 at 159
%     S  = transpose(window_1D_vec_around_center_col(transpose(S),round(length(S)/2),4.8,30,0,0));
    
	if (SR_zero_init_guess)
		f0 = zeros(N,1);                                                % N
	else
		f0 = transpose(interp1(1:M,abs(S),linspace(1,M,N)));            % N
	end;
	
	% (2.2)  Iteratively calculate the SR vector
	new_S = increase_1Dvector_res(f0,S,A,M,y_axis,SRnIter);

	if (DEBUG_FLAG >= 2) && plot_f
		figure;
		plot(linspace(y_axis(1),y_axis(end),M),abs(S)/max(abs(S))        ,'k.-'); hold on;
		plot(linspace(y_axis(end),y_axis(1),N),abs(new_S)/max(abs(new_S)),'r.-');
		legend({'Initial S','New S'});

		taPE_tmp = TaPE*linspace(0,length(y_axis)-1,length(y_axis));
		k_tmp    = gammaHz*GaPE*taPE_tmp;
		ph_tmp   = -(alpha2*(y_axis.^2) - alpha1*y_axis - alpha0) + k_tmp.*y_axis;
		figure; subplot(4,1,1); plot(abs(S)                    ,'b.-'); title('abs S');
		        subplot(4,1,2); plot(phase(S)                  ,'b.-'); title('phase S');
		        subplot(4,1,3); plot(abs(fftshift(fft(S)))     ,'b.-'); title('FFT S');
		        subplot(4,1,4); plot(abs(fftshift(fft(abs(S)))),'b.-'); title('FFT |S|');
		S = transpose(S) .* exp(-i*2*pi * ph_tmp);
		        subplot(4,1,1); hold on;plot(abs(S)                    ,'r-'); legend({'S','S*e^{-i*\phi}'});
		        subplot(4,1,2); hold on;plot(phase(S)                  ,'r-'); legend({'S','S*e^{-i*\phi}'});
		        subplot(4,1,3); hold on;plot(abs(fftshift(fft(S)))     ,'r-'); legend({'S','S*e^{-i*\phi}'});
		        subplot(4,1,4); hold on;plot(abs(fftshift(fft(abs(S)))),'r-'); legend({'S','S*e^{-i*\phi}'});
	end;

	% Filter out one peak
	if (SR_cs.exp)
		t = linspace(0,Ta,length(new_S));
		a1 = new_S;
		cen_col = round(N/2);
		Dk      = abs(NSE/Ly); % [1/cm]
		Dk_Hz   = Dk*Ly/Ta;    % [Hz]
		Dcs_Hz  = 0;           % [Hz]

		sft = SR_cs.sft;
		c2  = SR_cs.c2;

		c1  = 0.5*Dcs_Hz - sft;
		a2 = a1 .* transpose(exp(2*pi*i*c1*t));
		plot_filter = 0;
		a3 = fftshift(fft(a2));

		if (DEBUG_FLAG >= 1) && (RO_idx == check_idx)
			plot_filter = 1;
% 			figure;plot(smooth(abs(a3),2),'k.-');
% 			figure; % subplot(2,1,1); plot(abs(a3),'k.-');axis([0 70 0 100000]);
% 			figure; plot(abs(fftshift(fft(zero_pad(a2)))),'k-');axis([0 128 0 85000]);
		end;
% 		figure; axx = gca; % ax_spect
		a4 = window_1D_vec_around_center_col(transpose(a3),cen_col,SR_cs.exp,SR_cs.ESP,plot_filter,ax_spect); plot_filter=0;
		a5 = a4 .* exp(-2*pi*i*(c1+c2)*t/Ta);
		a6 = fftshift(fft(fftshift(a5)));
		new_S = a6(end:-1:1);
	else
		if (SR_cs.c2)
			FFT_new_S = fftshift(fft(new_S));
			FFT_new_S = FFT_new_S .* transpose(exp(2*pi*i*(SR_cs.c2)*(1:length(new_S))));
			new_S     = ifftshift(ifft(FFT_new_S));
		end;
	end;
	
	% (2.3)  Set the super resolution-ed solution into a SR matrix
	sig_mat_both_SR(RO_idx,:) = new_S;
end;

if (exist('pix_shift_ax','var') && (pix_shift_ax ~= 0) && SR_fix_inhomo_flag)
	axes(pix_shift_ax);
	imagesc(mapX,mapY,transpose(pixel_shift_mat)); set(gca,'ydir','normal'); colorbar; axis image;
end;

% if (SR_zero_init_guess)
sig_mat_both_SR = flipdim(sig_mat_both_SR,2); % don't know why this is needed but it seems that somtimes it is... :)
% end;

return;

