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
function [sig_mat_both_SR]=increase_2D_Hybrid_image_resolution4(sig_mat_both,Yie,Yfe,Ly,initial_dy,required_dy,...
                                              alpha0,alpha1,alpha2,GePE,Tp,NSE,TaPE,Ta,GaPE,SE_flag,SRExp,SRESP,...
                                              SR_ones_win_flag,SRnIter,SR_zero_init_guess,SR_fix_inhomo_flag,...
                                              SR_cs,ax_spect,T2PE_factor,pix_shift_ax)
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
Wcs       = -1e3;
CS_factor = 1;
te        = linspace(0,Tp,length(y_axis));
te_HR     = linspace(0,Tp,length(y_axis)*HRF);

% Excitation phase
a0 = 2*pi*alpha0;                                             % [rad]
a1 = 2*pi*alpha1;                                             % [rad/cm]
a2 = 2*pi*alpha2;                                             % [rad/cm^2]
phi_e    = (a2*(   y_axis.^2) + a1*   y_axis + a0);           % [rad]    N      Excitation phase axis
phi_e_HR = (a2*(HR_y_axis.^2) + a1*HR_y_axis + a0);           % [rad]    NxHRF  High resolution ...

phi_e_cs    = 2*pi*Wcs*te;
phi_e_cs_HR = 2*pi*Wcs*te_HR;

% if (SE_flag), GaPE  = -GaPE;                      end;      % (only for some multi-SP fid's: 43,...)
if (SE_flag), phi_e    = -phi_e;    phi_e_HR    = -phi_e_HR;    end;
if (SE_flag), phi_e_cs = -phi_e_cs; phi_e_cs_HR = -phi_e_cs_HR; end;

taPE = TaPE*linspace(0,M-1,M);                                % [sec]    Acquisition temporal axis
k    = 2*pi * gammaHz*GaPE*taPE;                              % [rad/cm]
k_cs = 2*pi*Wcs*taPE; % [rad]
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

A    = exp(i*(ones(M,1)*phi_e    + transpose(k)*y_axis));     % [rad]    MxN (M lines, N columns)
A_HR = exp(i*(ones(M,1)*phi_e_HR + transpose(k)*HR_y_axis));  % [rad]    Mx(N*HRF) (M lines, N*HRF columns)

A_cs    = exp(i*(ones(M,1)*phi_e_cs    + transpose(k_cs)*ones(1,N)));
A_cs_HR = exp(i*(ones(M,1)*phi_e_cs_HR + transpose(k_cs)*ones(1,N*HRF)));

% A    = A .* CS_factor .* A_cs;
% A_HR = A_HR .* CS_factor .* A_cs_HR;

% ------------------------
%  Load inhomogeneity map 
% ------------------------
if (SR_fix_inhomo_flag)
	HRF_B0  = HRF;
	img_LPE = 3.0;
	img_LRO = 2.8;
	map_LPE = 3.0;
	map_LRO = 2.8;
	map_fnm = 'XYmap_12may09_9.mat';
	[map_mat,mapX,mapY,te,ta] = SR_load_inhomo_map(map_fnm,sig_mat_both,Tp,Ta,M,N,...
	                                               img_LPE,img_LRO,map_LPE,map_LRO);

	% Calculate High-res map
	mapY_HR = linspace(mapY(1),mapY(end),length(mapY)*HRF);
	[xi,yi] = meshgrid(linspace(1,size(map_mat,2),size(map_mat,2)*HRF_B0), ...
	                   linspace(1,size(map_mat,1),size(map_mat,1)));
	map_mat_HR = interp2(map_mat,xi,yi);
	if (DEBUG_FLAG >= 2), figure; imagesc(transpose(map_mat_HR)); title('High-res map mat'); set(gca,'YDir','normal'); end; % caxis([-300 300]);

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
figure;
subplot(1,2,1); imagesc(ones(M,1)*phi_e_HR + transpose(k)*HR_y_axis); set(gca,'YDir','normal'); title('A [rad]');
subplot(1,2,2); imagesc(W_f);                                         set(gca,'YDir','normal'); title('W f [normalized]');
% figure;
% subplot(1,2,1); imagesc(ones(M,1)*phi_e_HR    + transpose(k)*HR_y_axis + ...
%                         ones(M,1)*phi_e_cs_HR + transpose(k_cs)*ones(1,N*HRF)); title('A [rad] w/ CS');
end;

if (~SR_fix_inhomo_flag)
	A = A .* W_f;
else
	A0 = A;
end;

% (2) Loop over all PE lines and apply the SR algorithm
check_idx = 20; %round(size(sig_mat_both,1)/2);
for RO_idx = 1:size(sig_mat_both,1)
	plot_f = 0;
	if (RO_idx == check_idx)
		disp(''); plot_f = 1;% check_idx = 90;
	end;

	if(SR_fix_inhomo_flag)
		% Extract the inhomogeneity vector
% 		dB0     = map_mat(RO_idx,end:-1:1);                    % [Hz] Probably the correct direction
% 		dB0_HR  = map_mat_HR(RO_idx,end:-1:1);                 %      Probably the correct direction
		dB0     = map_mat(RO_idx,:);                           % [Hz] Probably incorrect
		dB0_HR  = map_mat_HR(RO_idx,:);                        % [Hz] Probably incorrect

		% Smooth the inhomogeneity pattern
% 		poly_order = 5;
% 		y   = y_axis(end:-1:1);
% 		a   = polyfit(y(dB0~=0),dB0(dB0~=0),poly_order);
% 		dB0 = polyval(a,y);
% 		
% 		y_HR   = HR_y_axis(end:-1:1);
% 		a_HR   = polyfit(y_HR(dB0_HR~=0),dB0_HR(dB0_HR~=0),poly_order);
% 		dB0_HR = polyval(a_HR,y_HR);
		
% 		dB0     = transpose(smooth(dB0,5));
% 		dB0_HR  = transpose(smooth(dB0_HR,5*HRF));
		
		% Calculate the phase term, accumulated due to inhomogeneity (excitation; acquisition; total)
		if (SE_flag)  Phi_e_sign = -1; else Phi_e_sign = +1; end;
		dPh_B0e    = 2*pi * (transpose(te) * Phi_e_sign*dB0);             % [rad]
		dPh_B0e_HR = 2*pi * (transpose(te) * Phi_e_sign*dB0_HR);
		dPh_B0a    = 2*pi * (transpose(ta) * dB0);                        % [rad]
		dPh_B0a_HR = 2*pi * (transpose(ta) * dB0_HR);
		dPh_B0     = dPh_B0e    + dPh_B0a;                                % [rad]
		dPh_B0_HR  = dPh_B0e_HR + dPh_B0a_HR;
% 		dPh_B0     = ones(M,1) * dPh_B0;
% 		dPh_B0_HR  = ones(M,1) * dPh_B0_HR;

		% Calculate the transformation matrix, corresponding to the field inhomogeneity
		A_dB0 = exp(i*dPh_B0);

		% Calculate weighting, owing to inhomogeneity-phase-variation intra pixel dephasing
		Ph_A0    = ones(M,1)*phi_e    + transpose(k)*y_axis;
		Ph_A0_HR = ones(M,1)*phi_e_HR + transpose(k)*HR_y_axis;
		dB0_W_f  = SR_calc_dB0_W_f(exp(i*(dPh_B0_HR + Ph_A0_HR)),N,HRF_B0,SRExp,SRESP,SR_ones_win_flag);

% 		dB0_W_f_0 = SR_calc_dB0_W_f(exp(i*(dPh_B0_HR)),N,HRF_B0,0,0,0,dB0);
% 		dB0_W_f   = dB0_W_f_0 .* W_f;
		dB0_W_f   = dB0_W_f / max(max(dB0_W_f));

		% Set the final transformation matrix form
		A = A0 .* A_dB0 .* dB0_W_f;
		
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
			subplot(2,2,1); imagesc(dPh_B0_HR           ); title('dPh B0 HR');            colorbar; set(gca,'YDir','normal');
			subplot(2,2,2); imagesc(Ph_A0_HR            ); title('Ph A0 HR');             colorbar; set(gca,'YDir','normal');
			subplot(2,2,3); imagesc(dPh_B0_HR + Ph_A0_HR); title('dPh B0 HR + Ph A0 HR'); colorbar; set(gca,'YDir','normal');
			subplot(2,2,4); imagesc(dB0_W_f             ); title('dB0 W f');              colorbar; set(gca,'YDir','normal');

% 			figure;         imagesc(dB0_W_f_0);            title('dB0 W f 0');            colorbar; set(gca,'YDir','normal');
		end;

		% Calculate a pixel shift map
		PE_axis_hat  = y_axis + dB0/(GePE*gammaHz);
		pixel_length = abs(Ly)/length(y_axis);
		pixel_shift_mat(RO_idx,:) = (PE_axis_hat - y_axis)/pixel_length;
	end;
	
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
	end;

	% Filter out one peak
	if (SR_cs.exp)
		t = linspace(0,Ta,length(new_S));
		a1 = new_S;
		cen_col = round(N/2);
		Dk      = abs(NSE/Ly); % [1/cm]
		Dk_Hz   = Dk*Ly/Ta;    % [Hz]
		Dcs_Hz  = 1000;        % [Hz]

		sft = SR_cs.sft;
		c2  = SR_cs.c2;

		c1  = 0.5*Dcs_Hz - sft;
		a2 = a1 .* transpose(exp(2*pi*i*c1*t));
		plot_filter = 0;
		if (DEBUG_FLAG >= 0) && (RO_idx == check_idx)
			plot_filter = 1;
		end;
		a3 = fftshift(fft(a2));
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

if (exist('pix_shift_ax','var') && SR_fix_inhomo_flag)
	axes(pix_shift_ax);
	imagesc(mapX,mapY,transpose(pixel_shift_mat)); set(gca,'ydir','normal'); colorbar; axis image;
end;

% if (SR_zero_init_guess)
sig_mat_both_SR = flipdim(sig_mat_both_SR,2); % don't know why this is needed but it seems that it is...
% end;
return;

% ======================================================================================================= %
% ======================================================================================================= %

function [map_mat,mapX,mapY,te,ta] = SR_load_inhomo_map(map_fnm,sig_mat,Tp,Ta,M,N,img_LPE,img_LRO,map_LPE,map_LRO)
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
map_mat = transpose(tmp.map_mat);
mapX    = tmp.x_axis;
mapY    = tmp.y_axis;

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

[xi,yi] = meshgrid(linspace(1,size(map_mat,2),N), ...
				   linspace(1,size(map_mat,1),size(sig_mat,1)));
map_mat = interp2(map_mat,xi,yi);
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
