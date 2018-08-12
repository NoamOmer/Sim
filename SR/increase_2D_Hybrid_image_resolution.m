
% Increase the resolution of a 2D Hybrid SE/TE image.
% Image should be is oversampled (M > N) in order for the algorithm to converge.

% ================
% == Parameters ==
% ================
% sig_real      :  vector  : [none] :  Real      part of signal
% sig_imag      :  vector  : [none] :  Imaginary part of signal
% plot_flag     :  scalar  : [none] :  Enable / disable plots
% UD            :  struct  :
% UD.Lx         :  scalar  : [cm]   :  RO axis length
% UD.fid_offset :  scalar  : [none] :  Positive / negative echoes offset
% UD.Recon...   :  cell    : [none] :  Final RO axis, SE axis & image matrix
% UD.FID_axes   :  handle  : [none] :  Handle to axes on which the FID will be plotted
% initial_dy    :  scalar  : [cm]   :  Initial  SE axis pixel size
% required_dy   :  scalar  : [cm]   :  Required SE axis pixel size

function [base_offset, UD] = increase_2D_Hybrid_image_resolution(sig_real,sig_imag,plot_flag,UD,npRO,npPE,NPE,...
                                                           Yia,Yfa,Lx,Ly,alpha0,alpha1,alpha2,TaPE,GaPE,phi_0,...
                                                           shearing_factor,shearNegFactor,smoothF,winF,winType,...
                                                           initial_dy,required_dy,SE_flag,context)
set_globals;

% -----------------------------------------------
% Reshape FID signal into a 2D martix
% -----------------------------------------------
sig_comp = sig_real + i*sig_imag;
sig_mat  = reshape(sig_comp,2*(npRO + npPE),NPE);

% Adjust the PE axis to be from - to +.
% If acquisition was performed from Y-final to Y-initial (which is usually the case),
% flip the FID matrix. The resulting matrix rows' order (1:end) will match the excitation direction
if (Yfa < Yia)
	sig_mat = flipdim(sig_mat,2);
end;
figure;imagesc(abs(sig_mat)); title('Sig mat')

% -----------------------------------------------
% Create image axes
% -----------------------------------------------
if (Yfa < Yia)
	Yia_ = Yfa;
	Yfa_ = Yia;
end;
RO_axis = linspace(-Lx/2 , +Lx/2 , npRO );
PE_axis = linspace( Yia_ ,  Yfa_ , 2*NPE);

% Plot the FID matrix
% increase_2D_Hybrid_image_resolution_plot_FID_matrix(sig_mat,plot_flag)

% Matrix colomns: sample PE-axis : plotted along y-axis of imagesc
% Matrix lines  : sample RO-axis : plotted along x-axis of imagesc

% --------------------------------------------------------
%  First order fix:
%  Change the negative echo series by a fraction-of-a-pixel
%  factoring the spatial domain with exp(i*k0*RO_Axis)
% --------------------------------------------------------
% Calculate the offset between positive and negative echos
npRmpRO = 0;
[dummy,max_idx] = max(abs(sig_mat(round(npRmpRO + npRO/2),:)));
echo_PE_idx = max_idx;
pos_echo_start = npRmpRO + 1;
pos_echo_end   = npRmpRO + npRO;
neg_echo_start = npRmpRO + npRO + npRmpRO + npPE + npRmpRO + npRO;
neg_echo_end   = npRmpRO + npRO + npRmpRO + npPE + npRmpRO + 1;

[dummy,pos_echo_idx] = max(smooth(abs(sig_mat((pos_echo_start  : pos_echo_end),echo_PE_idx))));
[dummy,neg_echo_idx] = max(smooth(abs(sig_mat((neg_echo_start:-1:neg_echo_end),echo_PE_idx))));
base_offset = neg_echo_idx - pos_echo_idx;

% sig_mat2 = sig_mat;
% 
% for idx = -5:5
% sig_mat = sig_mat2;
% UD.fid_offset = idx;
if (UD.fid_offset)
	sig_mat = transpose(sig_mat);
	dkRO    = 1/Lx;	
	for(idx=1:NPE)
		fid_mat_neg(idx,:) = sig_mat(idx,neg_echo_start:-1:neg_echo_end);
	end;
	hybrid_imag_neg = fftshift(fft(fftshift(fid_mat_neg,2),[],2),2);
	hybrid_imag_neg = hybrid_imag_neg .* (ones(NPE,1) * exp(i*(2*pi*0.1*dkRO*UD.fid_offset)*RO_axis));  % moving in steps of 0.1 pixel
	fid_mat_neg     = ifftshift(ifft(ifftshift(hybrid_imag_neg,2),[],2),2);
	for (idx=1:NPE)
		sig_mat(idx,neg_echo_start:-1:neg_echo_end) = fid_mat_neg(idx,:);
	end;
	sig_mat = transpose(sig_mat);
	figure;imagesc(abs(sig_mat)); title('Sig mat - post offset adjustment')
end;

% -----------------------------------------------
% Combine positive and negative echoes
% -----------------------------------------------
for(idx=1:NPE)
	sig_mat_both(:,idx*2  ) = sig_mat((pos_echo_start:pos_echo_end),idx);
	sig_mat_both(:,idx*2-1) = sig_mat(((neg_echo_start-base_offset):-1:(neg_echo_end-base_offset)),idx);
end;
if (DEBUG_FLAG >= 0)
	figure;imagesc(abs(sig_mat_both)); title('Sig mat both - post combining')
	% increase_2D_Hybrid_image_resolution_plot_merged_FID_matrix(UD,sig_mat_both,plot_flag);
end;
% end;

% ----------------------------------------------------------------------------------------
% Fix frequency domain linear phase, caused by field inhomogeneity, resulting in shift
% in spatial domain. Shift alternating lines on the image on the RO axis. Do the shifting
% indirectly, i.e., in frequency space.
% ----------------------------------------------------------------------------------------
if (phi_0)
	shift = phi_0*linspace(0,1,size(sig_mat_both,2));
	for (idx=1:NPE)
		sig_mat_both(idx*2-1,:) = sig_mat_both(idx*2-1,:) .* exp(-i*shift);
		sig_mat_both(idx*2  ,:) = sig_mat_both(idx*2  ,:) .* exp(+i*shift);
	end;
	figure;imagesc(abs(sig_mat_both)); title('Sig mat both - \Phi_0 correction')
end;

% -------------------------------------------------
% Shear the FID matrix
% -------------------------------------------------
if (shearing_factor)
	tmp_PE_axis = linspace(Yia,Yfa,size(sig_mat_both,2));
	sig_mat_both = shear(sig_mat_both,tmp_PE_axis,shearing_factor,0);
	figure;imagesc(abs(sig_mat_both)); title('Sig mat both - post shearing')
end;

if (shearNegFactor)
	for (idx=1:NPE)
		sig_mat_neg(idx,:) = sig_mat_both(2*idx,:);
	end;
	tmp_PE_axis = linspace(Yia,Yfa,size(sig_mat_neg,2));
	sig_mat_neg = shear(sig_mat_neg,tmp_PE_axis,shearNegFactor,0);
	for(idx=1:NPE)
		sig_mat_both(2*idx,:) = sig_mat_neg(idx,:);
	end;
	figure;imagesc(abs(sig_mat_both)); title('Sig mat both - post neg shearing')
end;

% ---------------------------------------------------------
% Weight the FID according to distance from the center echo
% ---------------------------------------------------------
if (winF)
	smoothed_mat = smooth_mat(sig_mat_both, smoothF);
	[max_val_arr,max_row_arr] = max(smoothed_mat);
	[max_val,max_col] = max(max_val_arr);
	max_row = max_row_arr(max_col);

	for(idx=1:size(sig_mat_both,1))
		sig_mat_both(idx,:) = window_1D_vec_around_center_col(sig_mat_both(idx,:),max_col,winF,5,0);
	end;
	figure;imagesc(abs(sig_mat_both)); title('Sig mat both - post weighting')
end;

% -----------------------------------------------
% 2nd Stage: FT along the RO axis
% -----------------------------------------------
RO_fft_done = 0;
if (1)
	% Zero pad
	[sig_mat_both,n_row,n_col] = zero_pad(transpose(sig_mat_both));
	sig_mat_both = transpose(sig_mat_both);

	% FT
	N = size(sig_mat_both,1);
	sig_mat_both = fftshift(fft(fftshift(sig_mat_both,1),N,1),1);
	RO_fft_done = 1;

	if (DEBUG_FLAG >= 1)
		figure; imagesc(RO_axis,PE_axis,transpose(abs(sig_mat_both)));
		title('post FT image'); set(gca,'ydir','normal');
	end;
end;

% ===============================================
%  Increase image resolution
% ===============================================
% Calculate the initial y-axis, matching the old resolution
initial_y = linspace(Yia_,Yfa_,round(Ly/initial_dy));

% Adjust the required-resolution value so that the PE axis will contain a round number of pixels
new_dy = Ly / round(Ly/required_dy);
HRF  = 100;                                                   % High-resolution factor
y_axis    = linspace(Yia_,Yfa_,round(Ly/new_dy)  );
HR_y_axis = linspace(Yia_,Yfa_,length(y_axis)*HRF);           % High-resolution x-axis: HRF pts per pixel

disp(sprintf('New pixel (%3.2f[mm]) is smaller than initial pixel (%3.2f[mm]) by a factor of %3.1f',new_dy*10,initial_dy*10,initial_dy/new_dy));

% ---------------------------------------
%  Calculate the transformation matrix A
% ---------------------------------------
N = length(y_axis);                                           %          Number of pixels on the PE axis
M = 2*NPE;                                                    %          Number acquisition points on the PE axis

% Excitation phase
a0 = 2*pi*alpha0;                                             % [rad]
a1 = 2*pi*alpha1;                                             % [rad/cm]
a2 = 2*pi*alpha2;                                             % [rad/cm^2]
phi_e    = (a2*(   y_axis.^2) + a1*   y_axis + a0);           % [rad]    N      Excitation phase axis
phi_e_HR = (a2*(HR_y_axis.^2) + a1*HR_y_axis + a0);           % [rad]    NxHRF  High resolution ...
if (DEBUG_FLAG >= 4)  figure;  subplot(3,1,1); plot(phi_e,'.-'); title('\Phi_e');  end;

% Acquisition phase
if (SE_flag)
	phi_e    = -phi_e;
	phi_e_HR = -phi_e_HR;
end;

ta = TaPE*linspace(M,0,M);                                    % [sec]    Acquisition temporal axis
k  = 2*pi * gammaHz*GaPE*ta; % (!) -->                        % [rad/cm] <-- (!) NOTE THAT 'k' IS REVERSED (in time, meaning that its first values are large and correspond to long Ga and its last value is 0)

if (DEBUG_FLAG >= 4)  
subplot(3,1,2); plot(k,'.-'); title('k');
subplot(3,1,3); plot(phi_e + k(end)*y_axis,'.-'); title('\Phi_e + k(t_{end})y');
figure;
subplot(1,2,1); imagesc(ones(M,1)*phi_e);     title('\phi_e(y)');
subplot(1,2,2); imagesc(transpose(k)*y_axis); title('\phi_a(y,t)');
end;

A    = exp(i*(ones(M,1)*phi_e    + transpose(k)*y_axis));     % [rad]    MxN (M lines, N columns)
A_HR = exp(i*(ones(M,1)*phi_e_HR + transpose(k)*HR_y_axis));  % [rad]    Mx(N*HRF) (M lines, N*HRF columns)

% ----------------------------------------
%  Perform the high-resolution evaluation
% ----------------------------------------
% (1)  Start by calculating the weights around the stationary point, at each time-point
%      This is done only once, and applied to all PE lines
W_f = calculate_spatio_temporal_weights(A_HR,N,M,HRF,0,0,15,0);

if (DEBUG_FLAG >= 3)
figure;
subplot(1,2,1); imagesc(ones(M,1)*phi_e_HR + transpose(k)*HR_y_axis); title('A  ');
subplot(1,2,2); imagesc(W_f);                                         title('W f');
end;
A = A .* W_f;

% (2) Loop over all PE lines and apply the SR algorithm
check_idx = 33;
for RO_idx = 1:size(sig_mat_both,1)
	% (2.1)  Extract an initial guess
	% - This might be done by FT-ing S, filtering all spectral components that are above
	%   the current resolution (meanining higher than 1/initial_dx) and IFT-ing the result.
	S  = transpose(sig_mat_both(RO_idx,:));                   % M
	f0 = transpose(interp1(1:M,abs(S),linspace(1,M,N)));      % N

	if (DEBUG_FLAG >= 4) && (RO_idx == check_idx)
		figure; plot(linspace(y_axis(1),y_axis(end),M), abs(S)  / max(abs(S)) ,'k.-'); hold on;
		        plot(y_axis                           , abs(f0) / max(abs(f0)),'r.-'); legend({'Original S','S0'});
	end;

	% (2.2)  Iteratively calculate the SR vector
	n_iter = 10;
	new_S = increase_1Dvector_res(f0,S,A,M,y_axis,n_iter);
	
	if (DEBUG_FLAG >= 4) && (RO_idx == check_idx)
		figure;
		plot(abs(S)/max(abs(S)),'k.-'); hold on; plot(linspace(1,length(S),length(new_S)),abs(new_S)/max(abs(new_S)),'r.-');
		legend({'Initial S','New S'});
	end;

	% (2.3)  Set the super resolution-ed solution into a SR matrix
	sig_mat_both_SR(RO_idx,:) = new_S;
end;

% -------------------------------------------------------
% 2nd Stage: FT along the RO axis (Post filtering option)
% -------------------------------------------------------
if (~RO_fft_done)
	% Zero pad
	[sig_mat_both,n_row,n_col] = zero_pad(transpose(sig_mat_both));
	sig_mat_both = transpose(sig_mat_both);
	[sig_mat_both_SR,n_row,n_col] = zero_pad(transpose(sig_mat_both_SR));
	sig_mat_both_SR = transpose(sig_mat_both_SR);

	% FT
	N = size(sig_mat_both,1);
	imag_both = abs(fftshift(fft(fftshift(sig_mat_both,1),N,1),1));
	N = size(sig_mat_both_SR,1);
	imag_both_SR = abs(fftshift(fft(fftshift(sig_mat_both_SR,1),N,1),1));
else
	imag_both = abs(sig_mat_both);
	imag_both_SR = abs(sig_mat_both_SR);
end;

if (DEBUG_FLAG >= 0)
	figure;
	subplot(1,2,1);
	imagesc(RO_axis,PE_axis,abs(transpose(imag_both)));    set(gca,'Ydir','normal');
	title(sprintf('Orig image')); xlabel('RO axis'); ylabel('PE axis');

	subplot(1,2,2);
	imagesc(RO_axis,PE_axis,abs(transpose(imag_both_SR))); set(gca,'Ydir','normal');
	title(sprintf('SR image')); xlabel('RO axis'); ylabel('PE axis');
end;

% -----------------------------------------------
% Set image on UserData
% -----------------------------------------------
UD.ReconstructedImage{1} = RO_axis;
UD.ReconstructedImage{2} = PE_axis;
UD.ReconstructedImage{3} = imag_both_SR;

return;


function increase_2D_Hybrid_image_resolution_plot_FID_matrix(sig_mat,plot_flag)
if (plot_flag)
    figure;
    imagesc(abs(transpose(sig_mat))); set(gca,'Ydir','normal');
    title(sprintf('Raw FID - Amp.')); xlabel('RO axis'); ylabel('PE axis'); set_gca;

    figure;
    subplot(2,2,1);
    imagesc(transpose(unwrap(angle(sig_mat),[],1))); set(gca,'Ydir','normal');  % DIM=1 unwrapping along cols
    title(sprintf('Raw FID - Phase (unwrapped along RO-axis)')); xlabel('RO axis'); ylabel('PE axis');
    subplot(2,2,2);
    imagesc(transpose(unwrap(angle(sig_mat),[],2))); set(gca,'Ydir','normal');  % DIM=2 unwrapping along rows
    title(sprintf('Raw FID - Phase (unwrapped along PE-axis)')); xlabel('RO axis'); ylabel('PE axis');
    subplot(2,2,3);
    imagesc(transpose(angle(sig_mat))); set(gca,'Ydir','normal');
    title(sprintf('Raw FID - angle')); xlabel('RO axis'); ylabel('PE axis');
    subplot(2,2,4);
    imagesc(transpose(abs(sig_mat))); set(gca,'Ydir','normal');
    title(sprintf('Raw FID - Amp. ')); xlabel('RO axis'); ylabel('PE axis');
end;
return;

function increase_2D_Hybrid_image_resolution_plot_merged_FID_matrix(UD,sig_mat_both,plot_flag)
if (plot_flag)
	figure;
	subplot(2,2,1);
	imagesc(transpose(unwrap(angle(sig_mat_both),[],1))); set(gca,'Ydir','normal');
	title(sprintf('Merged FID - Phase (unwrapped along RO-axis)')); xlabel('RO axis'); ylabel('PE axis');
	subplot(2,2,2);
	imagesc(transpose(unwrap(angle(sig_mat_both),[],2))); set(gca,'Ydir','normal');
	title(sprintf('Merged FID - Phase (unwrapped along PE-axis)')); xlabel('RO axis'); ylabel('PE axis');
	subplot(2,2,3);
	imagesc(angle(transpose(sig_mat_both))); set(gca,'Ydir','normal');
	title(sprintf('Merged FID - angle')); xlabel('RO axis'); ylabel('PE axis');
	subplot(2,2,4);
	imagesc(abs(transpose(sig_mat_both))); set(gca,'Ydir','normal');
	title(sprintf('Merged FID - Amp. ')); xlabel('RO axis'); ylabel('PE axis');

	axes(UD.FID_axes);
	imagesc(abs(transpose(sig_mat_both))); set(gca,'Ydir','normal');
	title(sprintf('Merged FID - Amp. ')); xlabel('RO axis'); ylabel('PE axis');
end;
return;

