% Ver 2 : increase_2D_Hybrid_image_resolution2
%       : Remove all FID matrix manipulations (should be done before invoking the function)
function [sig_mat_both_SR] = increase_2D_Hybrid_image_resolution2(sig_mat_both,Yie,Yfe,Yia,Yfa,Ly,initial_dy,required_dy,...
                                                                  alpha0,alpha1,alpha2,NPE,TaPE,GaPE)
set_globals;

y_axis = linspace(Yie,Yfe,2*NPE);

% Adjust the PE axis to be from - to +.
% If acquisition was performed from Y-final to Y-initial (which is usually the case),
% flip the FID matrix. The resulting matrix rows' order (1:end) will match the excitation direction
if (Yfa < Yia)
	sig_mat_both = flipdim(sig_mat_both,2);
	Yia_ = Yfa;
	Yfa_ = Yia;
else
	Yia_ = Yia;
	Yfa_ = Yfa;
end;

% Calculate the initial y-axis, matching the old resolution
initial_y = linspace(Yia_,Yfa_,round(Ly/initial_dy));

% Adjust the required-resolution value so that the PE axis will contain a round number of pixels
new_dy = Ly / round(Ly/required_dy);
HRF  = 100;                                                   % High-resolution factor
y_axis    = linspace(Yia_,Yfa_,round(Ly/new_dy)  );
HR_y_axis = linspace(Yia_,Yfa_,length(y_axis)*HRF);           % High-resolution x-axis: HRF pts per pixel

% disp(sprintf('New pixel (%3.2f[mm]) is smaller than initial pixel (%3.2f[mm]) by a factor of %3.1f',new_dy*10,initial_dy*10,initial_dy/new_dy));

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
ta = TaPE*linspace(M,0,M);                                    % [sec]    Acquisition temporal axis
k  = 2*pi * gammaHz*GaPE*ta; % (!) -->                        % [rad/cm] <-- (!) NOTE THAT 'k' IS REVERSED (in time, meaning that its first values are large and correspond to long Ga and its last value is 0)

if (DEBUG_FLAG >= 4)  
subplot(3,1,2); plot(k,'.-'); title('k');
subplot(3,1,3); plot(phi_e + k(1)*y_axis,'.-'); title('\Phi_e + ky');
end;

A    = exp(i*(ones(M,1)*phi_e    + transpose(k)*y_axis));     % [rad]    MxN (M lines, N columns)
A_HR = exp(i*(ones(M,1)*phi_e_HR + transpose(k)*HR_y_axis));  % [rad]    Mx(N*HRF) (M lines, N*HRF columns)

% ----------------------------------------
%  Perform the high-resolution evaluation
% ----------------------------------------
% (1)  Start by calculating the weights around the stationary point, at each time-point
%      This is done only once, and applied to all PE lines
W_f = calculate_spatio_temporal_weights(A_HR,N,M,HRF,0);

if (DEBUG_FLAG >= 3)
figure;
subplot(1,2,1); imagesc(ones(M,1)*phi_e_HR + transpose(k)*HR_y_axis); title('A  ');
subplot(1,2,2); imagesc(W_f);                                         title('W f');
end;

A = A .* W_f;

% (2) Loop over all PE lines and apply the SR algorithm
check_idx = 68;
for RO_idx = 1:size(sig_mat_both,1)
	% (2.1)  Extract an initial guess
	% - This might be done by FT-ing S, filtering all spectral components that are above
	%   the current resolution (meanining higher than 1/initial_dx) and IFT-ing the result.
	S  = transpose(sig_mat_both(RO_idx,:));                   % M
	f0 = transpose(interp1(1:M,abs(S),linspace(1,M,N)));      % N

	if (DEBUG_FLAG >= 4) || (RO_idx == check_idx)
		figure; plot(linspace(y_axis(1),y_axis(end),M), abs(S)  / max(abs(S)) ,'k.-'); hold on;
		        plot(y_axis                           , abs(f0) / max(abs(f0)),'r.-'); legend({'Original S','S0'});
		axis([min(y_axis),max(y_axis),0,1]);
	end;

	% (2.2)  Iteratively calculate the SR vector
	n_iter = 5;
	new_S = increase_1Dvector_res(f0,S,A,M,y_axis,n_iter);
	
	if (RO_idx == check_idx) || (DEBUG_FLAG >= 4)
		figure;
		plot(abs(S)/max(abs(S)),'k.-'); hold on; plot(linspace(1,length(S),length(new_S)),abs(new_S)/max(abs(new_S)),'r.-');
		legend({'Initial S','New S'});
	end;

	% (2.3)  Set the super resolution-ed solution into a SR matrix
	sig_mat_both_SR(RO_idx,:) = new_S;
end;

return;

