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
    Yie,Yfe,Ly,initial_dy,required_dy,alpha0,alpha1,alpha2,GePE,Tp,NSE,TaPE,Ta,GaPE,...
    SE_flag,SRExp,SRESP,SR_ones_win_flag,SRnIter,SR_HRF,SR_zero_init_guess,...
	SR_fix_inhomo_flag,SR_cs,ax,T2PE_factor,...
	pix_shift_ax,SR_resp_win,skip_initial_proc,skip_initia_proc_fn) %#ok<INUSL>

set_globals;
cond_num = 0;
% if (SRnIter ~= 1)
% 	uiwait(msgbox('# of iterations ~= 1 .... Make sure that "increase_1Dvector_res" function is not configured to iter=1 mode !'));
% end;

if (~skip_initial_proc)

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

HRF  = SR_HRF;                                       % High-resolution factor
HR_y_axis = linspace(Yie,Yfe,length(y_axis)*HRF);    % High-resolution x-axis: HRF pts per pixel

% disp(sprintf('New pixel (%3.2f[mm]) is smaller than initial pixel (%3.2f[mm]) by a factor of %3.1f',new_dy*10,initial_dy*10,initial_dy/new_dy));

% ---------------------------------------
%  Calculate the transformation matrix A
% ---------------------------------------
% CODE FOR RASER / HYBRID processing
fudge_flip = 1;
if (SE_flag == 2) || fudge_flip
%     GaPE = -GaPE;
%     alpha2 = -alpha2;
%     alpha1 = -alpha1;
%     alpha0 = -alpha0;
end;

N = length(y_axis);                                           %          Number of pixels (required from the SR) on the PE axis
M = NSE;                                                      %          Number acquisition time-points on the PE axis

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

% taPE = TaPE*linspace(1,M,M);                                % [sec]    Acquisition temporal axis
taPE = TaPE*linspace(0,M-1,M);                                % [sec]    Acquisition temporal axis
k    = 2*pi * gammaHz*GaPE*taPE;                              % [rad/cm]
dk   = abs(k(2) - k(1));                                      % [rad/cm]
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

A    = exp(1i*(ones(M,1)*phi_e    + transpose(k)*y_axis));     % [rad]    MxN (M lines, N columns)
A_HR = exp(1i*(ones(M,1)*phi_e_HR + transpose(k)*HR_y_axis));  % [rad]    Mx(N*HRF) (M lines, N*HRF columns)

% A_cs    = exp(i*(ones(M,1)*phi_e_cs    + transpose(k_cs)*ones(1,N)));
% A_cs_HR = exp(i*(ones(M,1)*phi_e_cs_HR + transpose(k_cs)*ones(1,N*HRF)));

% A    = A .* CS_factor .* A_cs;
% A_HR = A_HR .* CS_factor .* A_cs_HR;


% -----------------------------------------------
%  Calculate the weights of the inversion matrix
% -----------------------------------------------
% (1)  Start by calculating the weights around the stationary point, at each time-point
%      This is done only once, and applied to all PE lines
% SR_factor = size(sig_mat_both,2)*initial_dy/Ly;          % [nbe oct 13 2011] commented
SR_factor = initial_dy / new_dy;                           % [nbe oct 13 2011] added

[W_f_PSF, W_f_gauss] = calculate_spatio_temporal_weights(A_HR,N,M,HRF,Ta,T2PE_factor,0,SRExp,SRESP,SR_ones_win_flag,SR_factor,dk,initial_dy,new_dy);

if (SR_resp_win == 1), W_f_plot = W_f_PSF; else W_f_plot = W_f_gauss; end;

if ((DEBUG_FLAG >= 4) || ax.ax_A)       , A_phase  = ones(M,1)*phi_e_HR + transpose(k)*HR_y_axis; end;
if ((DEBUG_FLAG >= 0) || ax.ax_A_Atag)  , B_PSF   = A.*W_f_PSF;                                   end;
if ((DEBUG_FLAG >= 0) || ax.ax_AW_AWtag), B_Gauss = A.*W_f_gauss;                                 end;
inv1 = abs((B_PSF')*B_PSF);      inv1 = inv1 / max(inv1(:));
inv2 = abs((B_Gauss')*B_Gauss);  inv2 = inv2 / max(inv2(:));
if (ax.ax_A),         axes(ax.ax_A);        imagesc(A_phase  ,'Parent',ax.ax_A       ); set(gca,'YDir','normal');       end;
if (ax.ax_W_f),       axes(ax.ax_W_f);      imagesc(W_f_plot ,'Parent',ax.ax_W_f     ); set(gca,'YDir','normal');       end;
if (ax.ax_A_Atag),    axes(ax.ax_A_Atag);   imagesc(inv1     ,'Parent',ax.ax_A_Atag  ); c=caxis; caxis(c/10); colorbar; end;
if (ax.ax_AW_AWtag),  axes(ax.ax_AW_AWtag); imagesc(inv2     ,'Parent',ax.ax_AW_AWtag); c=caxis; caxis(c/10); colorbar; end;
if (ax.fh_PSF),       PSF_recon_model_PSF   = abs((B_PSF'  )*(A.*W_f_PSF));
                      Gauss_recon_model_PSF = abs((B_Gauss')*(A.*W_f_PSF));
                      figure(ax.fh_PSF);  subplot(221); imagesc(PSF_recon_model_PSF  ); axis image;       title('PSF   model PSF'); 
                                          subplot(222); imagesc(Gauss_recon_model_PSF); axis image;       title('Guass model PSF');
                                          subplot(223); plot(PSF_recon_model_PSF  (round(NSE/2),:),'.-'); title('PSF   model 1D profile');
                                          subplot(224); plot(Gauss_recon_model_PSF(round(NSE/2),:),'.-'); title('Guass model 1D profile'); end;
% figure; plot(inv1(size(inv1,1)/2,:),'.-'); title(sprintf('PSF (Ones flag = %1.0f)',SR_ones_win_flag));
if (DEBUG_FLAG >= 4)
figure;
subplot(211); plot(inv1(size(inv1,1)/2,:),'.-'); title('PSF');
subplot(212); plot(inv2(size(inv2,1)/2,:),'.-'); title('Gauss');

figure;
subplot(2,2,1); imagesc(A_phase); set(gca,'YDir','normal');
subplot(2,2,2); imagesc(W_f_PSF); set(gca,'YDir','normal');
subplot(2,2,3); imagesc(inv1);
subplot(2,2,4); imagesc(inv2);
end;

% (2) Calculate the inversion matrix condition number
if (0) % skip this calculation for now... (for speeding purposes)
% fh_rank = 0;
% fh_rank = figure;
if (ax.ax_EIG), axes(ax.ax_EIG);delete(get(gca,'children')); hold on; end;
for idx = 1:2
	if (idx == 1)
		A_recon = B_PSF;
		str = 'PSF';
	else
		A_recon = B_Gauss;
		str = 'Gauss';
	end;
		
	OPTS.disp = 0;
% 	eig_max = eigs(A_recon,1,'LM',OPTS); % Op1
% 	eig_min = eigs(A_recon,1,'SM',OPTS); % Op1
% 	eig_max = eigs(A_recon,1,'LA',OPTS);
% 	eig_min = eigs(A_recon,1,'SA',OPTS);
% 	eig_vec = eig(A_recon);              % Op2 (identical to Op1)
	% figure;plot(eig_vals1,'.-');
% 	eig_vec_max = max(eig_vec);
% 	eig_vec_min = min(eig_vec);
	
	norm_A        = norm(B_PSF,2);         % the sampling matrix is pre-determined during the experiment and will always be B_PSF
	norm_invA     = norm(A_recon',2);      % the deconvolution matrix is determined according to the chosen reconstruction method
% 	norm_fro_A    = norm(A_recon,'fro');
% 	norm_fro_invA = norm(A_recon,'fro');
		
% 	cond_num1  = abs(eig_max     / eig_min);
% 	cond_num2  = abs(eig_vec_max / eig_vec_min);
% 	cond_num3  = cond(A_recon);  % which is:  norm(A) * norm(pinv(A)) % seems most reasonable - cond of the reconstruction matrix
	cond_num4  = norm_A     * norm_invA;                              % same principle as 'cond' but using our inverted matrix
% 	cond_num4  = norm_invA;
% 	cond_num5  = norm_fro_A * norm_fro_invA;
	
	cond_num(idx) = cond_num4;
% 	cond_num(idx) = rank(A_recon) / size(A_recon,1);
% 	rank(A_recon)
	
	% Estimate the maximal SR factor
	stl={'k.-','b.-'};
	if (ax.ax_EIG)
		[X]=svd(A_recon);
		plot(X,stl{idx});
% 		cond_num(idx) = max(X) / min(X); % similar to "cond_num3 = cond(A_recon)" above
		a = axis; axis([a(1) length(X)+50 0 max(X)+10]);
% 		set(gca,'xticklabel',[]);
% 		subplot(2,2,idx*2-1); imagesc(X /(max(max(X))));
% 		title(sprintf('0.01 x EIGs %s',str));
% 		colorbar; caxis([0 0.01]);
% 		subplot(2,2,idx*2);   plot(diag(X,0),'.-');  % axis([0 0.3*size(A_recon,2) 0 max(diag(X,0))])
% 		title(sprintf('EIGs %s (Rank = %1.1f)',str,rank(B_PSF)));  % the rank should be calculated for the sampling matrix only (?)
	end
	if (ax.ax_EIG), legend({'p','g'}); end;
end;
end;

% % % % % % % % % if (1)
% % % % % % % % % 	for idx = 1:M
% % % % % % % % % 	v = A_HR(idx,:);                                              % N*HRF
% % % % % % % % % 	m = reshape(v,HRF,N);                                         % NxHRF
% % % % % % % % % 	w_f = mean(m);                                                % N
% % % % % % % % % 
% % % % % % % % % 	SRExp = 1;
% % % % % % % % % 	SRESP = 2;
% % % % % % % % % 	if (SRExp)
% % % % % % % % % 		[max_val,max_idx] = max(w_f);
% % % % % % % % % 		w_f = window_1D_vec_around_center_col(w_f,max_idx,SRExp,SRESP,0);
% % % % % % % % % 	end;
% % % % % % % % % 
% % % % % % % % % 	W_f(idx,1:N) = w_f;                                           % MxN
% % % % % % % % % 	end;
% % % % % % % % % 	W_f = abs(W_f) / max(max(abs(W_f)));
% % % % % % % % % 	W_f = flipdim(W_f,2);
% % % % % % % % % end;

if (SR_resp_win == 1)
	W_f = W_f_PSF;
else
	W_f = W_f_gauss;
end;

if (DEBUG_FLAG >= 2)
	figure;
	subplot(2,2,1); imagesc(A_phase);    set(gca,'YDir','normal'); title('Phase of A [rad]');
	subplot(2,2,2); imagesc(W_f);        set(gca,'YDir','normal'); title('W f [normalized]');
	subplot(2,2,3); imagesc(abs((A')*A));
end;


% if (~SR_fix_inhomo_flag)
	A  = A .* W_f;
% else
% 	A0 = A .* W_f;
% end;

% Save info for later processing
if (~isempty(skip_initia_proc_fn))% && ~(skip_initial_proc == 0))
	save(skip_initia_proc_fn,'N','M','A','y_axis','cond_num','dk','new_dy')
end;

else  % skip_initia_proc (far far above)
	load(skip_initia_proc_fn);
end;

% ----------------------------------------------------
%  Loop over all RO points and apply the SR algorithm
% ----------------------------------------------------
% REG_mat  = eye(size(A,2));
% REG_mat  = eye(size(A,2)*2) + circshift(eye(size(A,2)*2),[0 1])*(-1);
% REG_mat  = REG_mat(65:64*3,65:64*3);
% % % REG_mat  = flipdim(REG_mat,2);
% REG_beta = 1e1;

check_idx = round(size(sig_mat_both,1)/2);
for RO_idx = 1:size(sig_mat_both,1)
	plot_f = 0;
	if (RO_idx == check_idx)
		disp(''); plot_f = 1;
	end;

	% (2.1)  Extract an initial guess
	% - This might be done by FT-ing S, filtering all spectral components that are above
	%   the current resolution (meanning higher than 1/initial_dx) and IFT-ing the result.
	S  = transpose(sig_mat_both(RO_idx,:));                             % M

	CS_filter_pre_SR = 0;
	CS_filter_done   = 0;
	if (CS_filter_pre_SR) && (SR_cs.exp || SR_cs.c2)
	S = transpose(increase_2D_Hybrid_image_resolution5_filter_1D_vector(SR_cs,S,length(S),NSE,Ly,Ta,ax.ax_spect,RO_idx,check_idx,DEBUG_FLAG,dk,initial_dy,new_dy));
	CS_filter_done = 1;
	end;
	
	if SR_zero_init_guess || 1
		f0 = zeros(N,1);                                                 % N
	else
		f0 = transpose(interp1(1:M,S,linspace(1,M,N)));                  % N
	end;
	
	% (2.2)  Iteratively calculate the SR vector
% 	new_S =         increase_1Dvector_res(f0, transpose(A)*S, transpose(A)*A + REG_beta*REG_mat, M,y_axis,SRnIter);
% 	new_S = flipdim(increase_1Dvector_res(f0, transpose(A)*S, transpose(A)*A + REG_beta*REG_mat, M,y_axis,SRnIter),1);
% 	new_S = increase_1Dvector_res(f0,S,A+REG_beta*REG_mat,M,y_axis,SRnIter);
	new_S = increase_1Dvector_res(f0,S,A,M,y_axis,SRnIter);

	if (DEBUG_FLAG >= 2) && plot_f
		figure;
		plot(linspace(y_axis(1),y_axis(end),M),abs(S)/max(abs(S))        ,'k.-'); hold on;
		plot(linspace(y_axis(end),y_axis(1),N),abs(new_S)/max(abs(new_S)),'r.-');
		legend({'Initial S','New S'});

		taPE_tmp = TaPE*linspace(0,length(y_axis)-1,length(y_axis));
		k_tmp    = gammaHz*GaPE*taPE_tmp;
		ph_tmp   = -(alpha2*(y_axis.^2) - alpha1*y_axis - alpha0) + k_tmp.*y_axis;
% 		ph_tmp   = -ph_tmp; % [nbe] just for testing purposes
		figure; subplot(4,1,1); plot(abs(S)                    ,'b.-'); title('abs S');
		        subplot(4,1,2); plot(phase(S)                  ,'b.-'); title('phase S');
		        subplot(4,1,3); plot(abs(fftshift(fft(S)))     ,'b.-'); title('FFT S');
		        subplot(4,1,4); plot(abs(fftshift(fft(abs(S)))),'b.-'); title('FFT |S|');
		S = transpose(S) .* exp(-1i*2*pi * ph_tmp);
		        subplot(4,1,1); hold on;plot(abs(S)                    ,'r-'); legend({'S','S*e^{-i*\phi}'});
		        subplot(4,1,2); hold on;plot(phase(S)                  ,'r-'); legend({'S','S*e^{-i*\phi}'});
		        subplot(4,1,3); hold on;plot(abs(fftshift(fft(S)))     ,'r-'); legend({'S','S*e^{-i*\phi}'});
		        subplot(4,1,4); hold on;plot(abs(fftshift(fft(abs(S)))),'r-'); legend({'S','S*e^{-i*\phi}'});
	end;

	% Filter out one peak - post SR
	if (~CS_filter_done) && (SR_cs.exp || SR_cs.c2)
	new_S = increase_2D_Hybrid_image_resolution5_filter_1D_vector(SR_cs,new_S,N,NSE,Ly,Ta,ax.ax_spect,RO_idx,check_idx,DEBUG_FLAG,dk,initial_dy,new_dy);
	new_S = new_S(end:-1:1);
	end;
	
	% (2.3)  Set the super resolution-ed vector into a SR matrix
	sig_mat_both_SR(RO_idx,:) = new_S;
end;

if (exist('pix_shift_ax','var') && (pix_shift_ax ~= 0) && SR_fix_inhomo_flag)
	axes(pix_shift_ax);
	imagesc(mapX,mapY,transpose(pixel_shift_mat)); set(gca,'ydir','normal'); colorbar; axis image;
end;

sig_mat_both_SR = flipdim(sig_mat_both_SR,2); % don't know why this is needed but it seems that somtimes it is... :)

return;

function [new_S] = increase_2D_Hybrid_image_resolution5_filter_1D_vector(SR_cs,new_S,N,NSE,Ly,Ta,ax_spect,RO_idx,check_idx,DEBUG_FLAG,dk,dz_LR,dz_HR)
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
	a2 = a1 .* transpose(exp(2*pi*1i*c1*t));
	plot_filter = 0;
	a3 = fftshift(fft(fftshift(a2)));
	
	if (DEBUG_FLAG >= 1) && (RO_idx == check_idx)
		plot_filter = 1;
		% figure;plot(smooth(abs(a3),2),'k.-');
		% figure; % subplot(2,1,1); plot(abs(a3),'k.-');axis([0 70 0 100000]);
		% figure; plot(abs(fftshift(fft(zero_pad(a2)))),'k-');axis([0 128 0 85000]);
	end;

	if (~isfield(SR_cs,'win')), SR_cs.win = 3; end;
	switch SR_cs.win
	case 1           % Gaussian window, similar to the one in "calculate weights" ans based on SPEN parameters
%		exp_arg = ( idx*dk*(dz_LR^2)  +  (-N:-1)*dz_HR ) / (dz_LR*2*pi*1); brought here from calculate weights...
		exp_arg = (linspace(-N/2,N/2,N) * dz_HR ) / (dz_LR*2*pi);
		exp_arg = exp_arg * SR_cs.exp;                   % this is not really needed but I added it just for sake of additional flexibility from GUI
		w_f_gauss = exp(-0.5 * (exp_arg.^SR_cs.ESP) );   % SR_cs.ESP should be = 2
		a4 = transpose(a3) .* w_f_gauss;

	case 2           % simple smooth filter
		h=msgbox('Using smooth filter OP'); pause(0.5); close(h);
		w_f_gauss = abs(create_smooth_filter(N,1,SR_cs.exp,10,SR_cs.ESP,0));
		a4 = transpose(a3) .* w_f_gauss;

	case 3           % orig
		a4 = window_1D_vec_around_center_col(transpose(a3),cen_col,SR_cs.exp,SR_cs.ESP,plot_filter,ax_spect);
	end;

	if (plot_filter)
		axes(ax_spect);
		delete(get(ax_spect,'Children'));
% 			figure;
		plot(abs(a3) / max(abs(a3)),'b.-'); hold on;
		plot(abs(a4) / max(abs(a4)),'g.-');
		if (exist('w_f_gauss','var') == 2), plot(abs(w_f_gauss) / max(abs(w_f_gauss)),'r'  ); end;
	end;
	plot_filter=0;
	
	a5 = a4 .* exp(-2*pi*1i*(c1+c2)*t/Ta);
	a6 = ifftshift(ifft(ifftshift(a5)));
	new_S = a6(end:-1:1) / length(new_S);
else
	if (SR_cs.c2)
		FFT_new_S = fftshift(fft(new_S));
		FFT_new_S = FFT_new_S .* transpose(exp(2*pi*1i*(SR_cs.c2)*(1:length(new_S))));
		new_S     = ifftshift(ifft(FFT_new_S));
	end;
end;
return;

% % For Thesis:
% A0_thesis = A .* W_f;
% B0_thesis = A .* W_f_gauss;
% figure; imagesc(abs(A0_thesis));              title('A');
% figure; imagesc(abs(B0_thesis));              title('B');
% figure; imagesc(abs((A0_thesis')*A0_thesis)); title('A^H x A');
% figure; imagesc(abs((B0_thesis')*A0_thesis)); title('B^H x A  (B=Gaussian)');
% 
