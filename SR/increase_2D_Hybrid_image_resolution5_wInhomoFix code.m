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
if (SE_flag), phi_e    = -phi_e;    phi_e_HR    = -phi_e_HR;    end;
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

% ------------------------
%  Load inhomogeneity map 
% ------------------------
if (SR_fix_inhomo_flag)
	HRF_B0  = 1; %HRF;
	img_LPE = 3.0;
	img_LRO = 2.8;
	map_LPE = 3.0;
	map_LRO = 2.8;
	map_fnm = 'XYmap_13may09_6.mat';
	[map_mat,map_mask,mapX,mapY,te,ta] = SR_load_inhomo_map(map_fnm,sig_mat_both,Tp,Ta,M,N,...
	                                                        img_LPE,img_LRO,map_LPE,map_LRO);

	% Calculate High-res map
	[xi,yi] = meshgrid(linspace(1,size(map_mat,2),size(map_mat,2)*HRF_B0), ...
	                   linspace(1,size(map_mat,1),size(map_mat,1)*HRF_B0));
	map_mat_HR = interp2(map_mat,xi,yi);

% 	poly_order = 6;
% 	y   = y_axis(end:-1:1);
% 	for RO_idx = 1:size(map_mat,1)
% 		dB0 = map_mat(RO_idx,:);
% 		a   = polyfit(y(dB0~=0),dB0(dB0~=0),poly_order);
% 		if (length(dB0(dB0~=0)) < 20)
% 			new_dB0 = map_mat(RO_idx,:);
% 		else
% 			new_dB0 = polyval(a,y);
% 		end;
% 		new_dB0(new_dB0 > max(dB0)) = max(dB0);
% 		new_dB0(new_dB0 < min(dB0)) = min(dB0);
% 		map_mat1(RO_idx,:) = new_dB0;
% 	end;
% 	figure; subplot(2,1,1); imagesc(mapX,mapY,transpose(map_mat));  title('pre  interp'); set(gca,'YDir','normal');xlabel('RO-axis [cm]');ylabel('PE-axis [cm]');axis image; colorbar; ca=caxis;
% 	        subplot(2,1,2); imagesc(mapX,mapY,transpose(map_mat1)); title('post interp'); set(gca,'YDir','normal');xlabel('RO-axis [cm]');ylabel('PE-axis [cm]');axis image; colorbar; % caxis(ca);
	
	dB0_std_mat_caxis = 8;
	if (1)
% 		if (SE_flag)  Phi_e_sign = -1; else Phi_e_sign = +1; end;
% 		t_of_y = (Tp - te)*Phi_e_sign + (Ta - ta);                         % [sec]
% 		t_of_y = t_of_y(end:-1:1);
		Fx = 0.5;
		Fy = 4;
		for row = 1:size(map_mat,1)
			RO_low_idx = max(1                 ,round(1 + HRF_B0*(row-1) - Fx*HRF_B0/2));
			RO_hig_idx = min(size(map_mat_HR,1),round(1 + HRF_B0*(row-1) + Fx*HRF_B0/2));
			for col = 1:size(map_mat,2)
				PE_low_idx = max(1                 ,round(1 + HRF_B0*(col-1) - Fy*HRF_B0/2));
				PE_hig_idx = min(size(map_mat_HR,2),round(1 + HRF_B0*(col-1) + Fy*HRF_B0/2));
				voxel = map_mat_HR(RO_low_idx:RO_hig_idx,PE_low_idx:PE_hig_idx);
				dB0_std_mat(row,col) = std(voxel(:));

%				dB0_std_mat(row,col) = var(voxel(:));

% 				t_acq = t_of_y(col);
% 				att(row,col) = abs(mean(exp(i*2*pi*voxel(:)*t_acq)));
			end;
		end;
		dB0_std_mat(dB0_std_mat >= dB0_std_mat_caxis) = 0;
		att = exp(-dB0_std_mat);% .* (ones(size(dB0_std_mat,1),1)*(linspace(0,1,size(dB0_std_mat,2)))));
		
% 		att(att < 0.7) = 1;
	else
		tmp = load('13may09_2_T2star_map_n_params');
		dB0_std_mat = tmp.T2mat1;
		dB0_std_mat = transpose(smooth_mat(dB0_std_mat,9));

		[xi,yi] = meshgrid(linspace(1,size(dB0_std_mat,2),size(map_mat,2)), ...
						   linspace(1,size(dB0_std_mat,1),size(map_mat,1)));
		dB0_std_mat = interp2(dB0_std_mat,xi,yi);

		figure;
		subplot(2,2,1); imagesc(mapX,mapY,transpose(dB0_std_mat));set(gca,'YDir','normal'); axis image; colorbar; title('T2* [ms]');
		dB0_std_mat = dB0_std_mat * 1e-3;     % [sec]
		dB0_std_mat = 1./dB0_std_mat;         % [Hz]
		subplot(2,2,2); imagesc(mapX,mapY,transpose(dB0_std_mat)); set(gca,'YDir','normal'); axis image; colorbar; title('T2* [Hz]');
		caxis([0 50]);
		dB0_std_mat(dB0_std_mat < 3) = max(max(dB0_std_mat)); % [ms]
	% 	dB0_std_mat = dB0_std_mat.*map_mask;
		subplot(2,2,3); imagesc(mapX,mapY,transpose(dB0_std_mat)); set(gca,'YDir','normal'); axis image; colorbar; title('T2* [Hz] masked');
		caxis([0 dB0_std_mat_caxis]);

		att = exp(-dB0_std_mat);% .* (ones(size(dB0_std_mat,1),1)*(linspace(0,1,size(dB0_std_mat,2)))));
	end;
	
	if (DEBUG_FLAG >= 2), figure;
		subplot(2,2,1); imagesc(mapX,mapY,transpose(map_mat));     title('Map mat');          set(gca,'YDir','normal'); axis image; colorbar;
		subplot(2,2,2); imagesc(mapX,mapY,transpose(map_mat_HR));  title('High-res map mat'); set(gca,'YDir','normal'); axis image; colorbar;
% 		caxis([]);
		subplot(2,2,3); imagesc(mapX,mapY,transpose(dB0_std_mat)); title('dB_0 STD mat');     set(gca,'YDir','normal'); axis image; colorbar;
		caxis([0 dB0_std_mat_caxis])
		subplot(2,2,4); imagesc(mapX,mapY,transpose(att));         title('Att (w/o temporal weighting)');set(gca,'YDir','normal'); axis image; colorbar;
% 		caxis([0 1])
	end;
% return;
	te_HR = linspace(te(1),te(end),length(te)*HRF_B0);
	ta_HR = linspace(ta(1),ta(end),length(ta)*HRF_B0);
end;

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

% (2) Loop over all PE lines and apply the SR algorithm
check_idx = round(size(sig_mat_both,1)/2); % 75;
for RO_idx = 1:size(sig_mat_both,1)
	plot_f = 0;
	if (RO_idx == check_idx)
		disp(''); plot_f = 1;
	end;

	if(SR_fix_inhomo_flag)
		% Extract the inhomogeneity vector
		dB0 = map_mat(RO_idx,:);                           % [Hz] Probably incorrect

		% Smooth the inhomogeneity pattern
% figure; plot(dB0,'.-'); hold on;
% 		poly_order = 2;
% 		y   = y_axis(end:-1:1);
% 		a   = polyfit(y(dB0~=0),dB0(dB0~=0),poly_order);
% 		dB0 = polyval(a,y);
% plot(dB0,'k.-');

% 		dB0     = transpose(smooth(dB0,5));
		
		% Calculate the phase term, accumulated due to inhomogeneity (excitation; acquisition; total)
		if (SE_flag)  Phi_e_sign = -1; else Phi_e_sign = +1; end;
		dPh_B0e    = 2*pi * (transpose(te) * Phi_e_sign*dB0);             % [rad]
		dPh_B0a    = 2*pi * (transpose(ta) * dB0);                        % [rad]
		dPh_B0     = dPh_B0e    + dPh_B0a;                                % [rad]

		% Calculate the transformation matrix, corresponding to the field inhomogeneity
		A_dB0 = exp(i*dPh_B0);

		% Calculate weighting, owing to inhomogeneity-phase-variation intra pixel dephasing
% 		if (RO_idx > 1)                   ,  dB0_pre  = map_mat(RO_idx-1,:);  else  dB0_pre  = dB0;  end;
% 		if (RO_idx < size(sig_mat_both,1)),  dB0_post = map_mat(RO_idx+1,:);  else  dB0_post = dB0;  end;
% 		df_x = ((dB0_post - dB0) + (dB0 - dB0_pre)) / 2;
% 		
% 		dB0_extend = [dB0(1) dB0 dB0(end)];
% 		df_y = ((dB0_extend(3:end) - dB0_extend(2:end-1)) + (dB0_extend(2:end-1) - dB0_extend(1:end-2)) / 2);
% 		
% 		df = sqrt(df_x.^2 + df_y.^2);
% 		t_of_y = (Tp - te)*Phi_e_sign + (Ta - ta);
% 		t_of_y = t_of_y(end:-1:1);
% 		dB0_W_f = exp(-transpose(t_of_y)*df);
% % 		dB0_W_f   = dB0_W_f / max(max(dB0_W_f));
% 
% dfx_mat(RO_idx,:) = df_x;
% dfy_mat(RO_idx,:) = df_y;
% df_mat(RO_idx,:)  = df;
% Dt_mat(RO_idx,:)  = t_of_y;
% att_mat(RO_idx,:) = exp(-df.*t_of_y);

		% Calculate weighting, owing to inhomogeneity-phase-variation intra pixel dephasing
% 		RO_low_idx = max(1                   ,round(RO_idx - HRF_B0/2));
% 		RO_hig_idx = min(size(sig_mat_both,1),round(RO_idx + HRF_B0/2));
% 		dB0_HR_voxel = map_mat_HR(RO_low_idx:RO_hig_idx,:);
% 		dB0_std = std(dB0_HR_voxel);                                       % [Hz]
		dB0_std_vec = dB0_std_mat(RO_idx,:);
		t_of_y = (Tp - te)*Phi_e_sign + (Ta - ta);                         % [sec]
		t_of_y = t_of_y(end:-1:1);
		dB0_W_f = exp(-transpose(t_of_y)*dB0_std_vec);

% dB0_std_mat(RO_idx) = dB0_std;
Dt_mat(RO_idx,:)    = t_of_y;
att_mat(RO_idx,:)   = exp(-dB0_std_vec.*t_of_y);% .* exp(-(t_of_y + 5e-3) / (500e-3));

% att_mat(RO_idx,:)   = att(RO_idx,:);

% T2_star = 20e-3;
% att_mat(RO_idx,:)   = exp(-t_of_y/T2_star);% .* exp(-t_of_y);

		% Set the final transformation matrix form
		A = A0;% .* A_dB0 .* dB0_W_f;
		
		if plot_f && (DEBUG_FLAG >= 2)
			figure;
			subplot(2,2,1); imagesc(transpose(abs(sig_mat_both)));             title('Sig');    set(gca,'ydir','normal');
			subplot(2,2,2); imagesc(transpose(map_mat));                       title('Map');    set(gca,'ydir','normal');
			subplot(2,2,3); plot(abs(transpose(sig_mat_both(RO_idx,:))),'.-'); title('S');
			subplot(2,2,4); plot(dB0,'.-');                                    title('dB0');

			figure;
			subplot(2,2,1); plot(dB0,'.-');                title('dB0    ');
			subplot(2,2,2); imagesc(dPh_B0e);              title('dPh B0e');              colorbar; set(gca,'YDir','normal');
			subplot(2,2,3); imagesc(dPh_B0a);              title('dPh B0a');              colorbar; set(gca,'YDir','normal');
			subplot(2,2,4); imagesc(dPh_B0);               title('dPh B0e + dPh B0a');    colorbar; set(gca,'YDir','normal');

			figure;
			subplot(2,2,1); plot(df_x,'.-');               title('df_x');
			subplot(2,2,2); plot(df_y,'.-');               title('df_y');
			subplot(2,2,3); plot(df,'.-');                 title('df');
			subplot(2,2,4); plot(y_axis,t_of_y,'.-');      title('\Deltat(y)');
			
			figure;
			subplot(2,2,1); imagesc(W_f);                  title('W f');            xlabel('y-axis'); ylabel('t_a'); colorbar; set(gca,'YDir','normal');
			subplot(2,2,2); imagesc(dB0_W_f);              title('dB0 W f');        xlabel('y-axis'); ylabel('t_a'); colorbar; set(gca,'YDir','normal');
			subplot(2,2,3); imagesc(dB0_W_f.*W_f);         title('dB0 W f .* W f'); xlabel('y-axis'); ylabel('t_a'); colorbar; set(gca,'YDir','normal');
			subplot(2,2,4); imagesc(angle(A));             title('angle(A)');       xlabel('y-axis'); ylabel('t_a'); colorbar; set(gca,'YDir','normal');
			
% 			figure;         imagesc(dB0_W_f_0);            title('dB0 W f 0');            colorbar; set(gca,'YDir','normal');
		end;
		
		% Calculate a pixel shift map
		PE_axis_hat  = y_axis + dB0/(GePE*gammaHz);
		pixel_length = abs(Ly)/length(y_axis);
		pixel_shift_mat(RO_idx,:) = (PE_axis_hat - y_axis)/pixel_length;
	end; % closing: if(SR_fix_inhomo_flag)
	
	% (2.1)  Extract an initial guess
	% - This might be done by FT-ing S, filtering all spectral components that are above
	%   the current resolution (meanining higher than 1/initial_dx) and IFT-ing the result.
	S  = transpose(sig_mat_both(RO_idx,:));                             % M
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
sig_mat_both_SR = flipdim(sig_mat_both_SR,2); % don't know why this is needed but it seems that it is...
% end;

if SR_fix_inhomo_flag
% dfx_mat(abs(dfx_mat)  > 10  ) = 0;
% dfy_mat(abs(dfy_mat)  > 10  ) = 0;
% df_mat (df_mat        > 15  ) = 0;
% % att_mat(att_mat       < 0.20) = 1;
% figure;
% subplot(2,3,1); imagesc(mapX,mapY,transpose(df_mat));               title('df MAT');  set(gca,'ydir','normal'); colorbar; axis image;
% subplot(2,3,2); imagesc(mapX,mapY,transpose(dfx_mat));              title('dfx MAT'); set(gca,'ydir','normal'); colorbar; axis image;
% subplot(2,3,3); imagesc(mapX,mapY,transpose(dfy_mat));              title('dfy MAT'); set(gca,'ydir','normal'); colorbar; axis image;
% % subplot(2,3,4); imagesc(mapX,mapY,transpose(Dt_mat));               title('Dt MAT');  set(gca,'ydir','normal'); colorbar; axis image; 
% subplot(2,3,4); imagesc(mapX,mapY,transpose(att_mat));              title('att MAT'); set(gca,'ydir','normal'); colorbar; axis image;
% subplot(2,3,5); imagesc(mapX,mapY,abs(transpose(sig_mat_both_SR))); title('Sig');     set(gca,'ydir','normal'); colorbar; axis image; cax = caxis;
% subplot(2,3,6); imagesc(mapX,mapY,flipdim(transpose(abs(sig_mat_both_SR ./ att_mat)),1)); caxis(cax); title('./ att'); colorbar; axis image; %

figure;
subplot(2,3,1); imagesc(mapX,mapY,transpose(dB0_std_mat));          title('dB0 std MAT');  set(gca,'ydir','normal'); colorbar; axis image;
caxis([0 dB0_std_mat_caxis]);
subplot(2,3,2); imagesc(mapX,mapY,transpose(Dt_mat));               title('Dt MAT');       set(gca,'ydir','normal'); colorbar; axis image;
subplot(2,3,3); imagesc(mapX,mapY,transpose(att_mat));              title('att MAT');      set(gca,'ydir','normal'); colorbar; axis image;
caxis([0.8 1]);
subplot(2,3,4); imagesc(mapX,mapY,transpose(map_mat_HR));           title('High-res map'); set(gca,'YDir','normal'); colorbar; axis image;
subplot(2,3,5); imagesc(mapX,mapY,abs(transpose(sig_mat_both_SR))); title('Sig');          set(gca,'ydir','normal'); colorbar; axis image;
cax = caxis;
subplot(2,3,6); imagesc(mapX,mapY,transpose(abs(sig_mat_both_SR./att_mat))); title('./att');set(gca,'ydir','normal'); colorbar; axis image;
caxis([0 3000]);
% caxis([0 50]);
end;
return;

% ======================================================================================================= %
% ======================================================================================================= %

function [map_mat,map_mask,mapX,mapY,te,ta] = SR_load_inhomo_map(map_fnm,sig_mat,Tp,Ta,M,N,img_LPE,img_LRO,map_LPE,map_LRO)
set_globals;

te = linspace(0,Tp,N);
ta = linspace(0,Ta,N);
imgX = linspace(-img_LRO/2,+img_LRO,size(sig_mat,1));
imgY = linspace(-img_LPE/2,+img_LPE,size(sig_mat,2));

% tmp = load(map_fnm); tmp.map_mat = tmp.mapHz; tmp.map_mask = ones(size(tmp.map_mat));
% load(map_fnm); %XYmap_6apr09_1; %4apr09_4_MAT %load('12may09_3_MAT');
% Shape3D = flipdim(Shape3D,2);
% tmp.map_mat  = Shape3D*1000;
% tmp.map_mask = ones(size(tmp.map_mat));
% tmp.x_axis   = Y;%(end:-1:1);
% tmp.y_axis   = X(end:-1:1);

% map_mat = tmp.map_mat .* tmp.map_mask;
% mapX    = tmp.x_axis;
% mapY    = tmp.y_axis;

tmp = load(map_fnm);
map_mat  = transpose(tmp.map_mat);
map_mask = transpose(tmp.map_mask);
mapX     = tmp.x_axis;
mapY     = tmp.y_axis;

% % Flip the y-axis to match the acquisition direction (from last y to first y)
% mapY    = flipdim(mapY,2);
% map_mat = flipdim(map_mat,2);

% if (DEBUG_FLAG >= 2)
% 	figure; subplot(1,2,1); imagesc(mapX,mapY,map_mat); title('Map w/  axes'); set(gca,'YDir','normal');
% 			subplot(1,2,2); imagesc(          map_mat); title('Map w/o axes'); set(gca,'YDir','normal'); end;

if (img_LPE ~= map_LPE) || (img_LRO ~= map_LRO)
	[map_ln,map_cl] = size(map_mat);

	m = zeros(round(map_ln*img_LRO/map_LRO),round(map_cl*img_LPE/map_LPE));
	[new_ln,new_cl] = size(m);
	low_row	 = abs(ceil((map_ln - new_ln)/2)) + 1;
	high_row = low_row + new_ln - 1;
	low_col	 = ceil((map_cl - new_cl)/2);
	high_col = low_col + new_cl - 1;
	sx = 0;
	sy = 0;
	m = map_mat(low_row+sx:high_row+sx,low_col+sy:high_col+sy);
	map_mat = flipdim(m,2);
	mapX = linspace(-img_LRO/2,+img_LRO/2,new_cl);
	mapY = linspace(-img_LPE/2,+img_LPE/2,new_ln);
end;

if (DEBUG_FLAG >= 2)
figure;
subplot(2,2,1);imagesc(imgX,imgY,transpose(abs(sig_mat)));title('Image w/  axes');set(gca,'YDir','normal');xlabel('RO-axis [cm]');ylabel('PE-axis [cm]');axis image;
subplot(2,2,2);imagesc(mapX,mapY,transpose(map_mat));     title('Map   w/  axes');set(gca,'YDir','normal');xlabel('RO-axis [cm]');ylabel('PE-axis [cm]');axis image;%caxis([-300 300]);
subplot(2,2,3);imagesc(          transpose(abs(sig_mat)));title('Image w/o axes');set(gca,'YDir','normal');xlabel('RO-axis [cm]');ylabel('PE-axis [cm]');
subplot(2,2,4);imagesc(          transpose(map_mat));     title('Map   w/o axes');set(gca,'YDir','normal');xlabel('RO-axis [cm]');ylabel('PE-axis [cm]');%caxis([-300 300]);
end;

% Interpolate / extrapolate uncharted map regions
	poly_order = 8;
	y   = 1:size(map_mat,1);
	for RO_idx = 1:size(map_mat,1)
		dB0 = map_mat(RO_idx,:);
		a   = polyfit(y(dB0~=0),dB0(dB0~=0),poly_order);
		if (length(dB0(dB0~=0)) < 10)
			new_dB0 = map_mat(RO_idx,:);
		else
			new_dB0 = polyval(a,y);
		end;
		new_dB0(new_dB0 > max(dB0)) = max(dB0);
		new_dB0(new_dB0 < min(dB0)) = min(dB0);
		map_mat1(RO_idx,:) = new_dB0;
	end;
	figure; subplot(2,1,1); imagesc(mapX,mapY,transpose(map_mat));  title('pre  interp'); set(gca,'YDir','normal');xlabel('RO-axis [cm]');ylabel('PE-axis [cm]');axis image; colorbar; ca=caxis;
	        subplot(2,1,2); imagesc(mapX,mapY,transpose(map_mat1)); title('post interp'); set(gca,'YDir','normal');xlabel('RO-axis [cm]');ylabel('PE-axis [cm]');axis image; colorbar; % caxis(ca);
% map_mat1 = map_mat;

% Interpolate the map to match the image resolution
[xi,yi] = meshgrid(linspace(1,size(map_mat1,2),N), ...
				   linspace(1,size(map_mat1,1),size(sig_mat,1)));
map_mat  = interp2(map_mat1,xi,yi);
map_mask = interp2(map_mask,xi,yi);
if (DEBUG_FLAG >= 3)
% subplot(2,2,3);
figure;imagesc(transpose(map_mat)); title('Low-res map mat'); set(gca,'YDir','normal'); % caxis([-300 300]);
end;
return;

% ======================================================================================================= %
% ======================================================================================================= %

function [dB0_W_f] = SR_calc_dB0_W_f(dB0_HR,N,HRF,SRExp,SRESP,SR_ones_win_flag,dB0)
for idx = 1:size(dB0_HR,1)
	v = dB0_HR(idx,:);
	v = reshape(v,HRF,N);
	w_f = abs(mean(v));                     % Calculate the intra pixel dephasing
% 	w_f = pi - abs(std(angle(v)));           % Calculate the intra pixel dephasing
	
	[max_val,max_idx] = max(w_f);           % find the current maximal weight
	if (SR_ones_win_flag)
		w_f = ones(1,N);                    % give identical weight to all points - this must be accompanied with SRExp~=0
	end;
	if (SRExp == 0 && SR_ones_win_flag) uiwait(msgbox('Warning: Using equal weight matrix w/o windowing')); error('Exiting'); end;

	% Limit the weighting response function around the SP region
	if (SRExp)
	w_f = window_1D_vec_around_center_col(w_f,max_idx,SRExp,SRESP,0);     % window around the incremented index
	end;

% 	w_f(dB0 == 0) = 1;
	
	W_f(idx,1:N) = w_f;
end;

dB0_W_f = abs(W_f) / max(max(abs(W_f)));    % Normalize the weighting function to [0..1]

return;
