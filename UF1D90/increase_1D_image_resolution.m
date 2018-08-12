
% Increase the resolution of a 1D image.
%
% S           :  1D vector  :  [none]     :  Measured signal     containint M elements
% a0          :     scalar  :  [rad]      :  Excitation phase constant  term coefficient
% a1          :     scalar  :  [rad/cm]   :  Excitation phase linear    term coefficient
% a2          :     scalar  :  [rad/cm^2] :  Excitation phase quadratic term coefficient
% Ga          :     scalar  :  [G/cm]     :  Acquisition gradient
% Ta          :     scalar  :  [sec]      :  Total acquisition time
% L           :     scalar  :  [cm]       :  Sample length
% initial_dx  :     scalar  :  [cm]       :  Initial  pixel size
% required_dx :     scalar  :  [cm]       :  Required pixel size
%
% Assumptions
% 1. Image is oversampled (M > N)

function [new_S] = increase_1D_image_resolution(S,a0,a1,a2,Ga,Ta,L,initial_dx,required_dx,SRExp,SRESP,SR_HRF)
set_globals;
% close all;

% Calculate the initial x-axis, matching the old resolution
initial_x = linspace(-L/2,+L/2,round(L/initial_dx));

% Adjust the required resolution so that it the sample length will contain a round number of pixels
new_dx = L / round(L/required_dx);
x = linspace(-L/2,+L/2,round(L/new_dx));

x_HR = linspace(-L/2,+L/2,length(x)*SR_HRF);          % High-resolution x-axis: HRF pts per pixel

disp(sprintf('New pixel (%3.1f[mm]) is smaller than initial pixel (%3.1f[mm]) by a factor of %3.1f',new_dx*10,initial_dx*10,initial_dx/new_dx));

% ---------------------------------------
%  Calculate the transformation matrix A
% ---------------------------------------
N = length(x);
M = length(S);

phi_e    = (a2*(x.^2)    + a1*x    + a0);                % [rad]    N      Excitation phase axis
phi_e_HR = (a2*(x_HR.^2) + a1*x_HR + a0);                % [rad]    NxHRF  High resolution ...

ta = linspace(0,Ta,M);                                   % [sec]    Note: ta goes from 0 to Ta and not reversed probably because S is also not flipped...
k  = -2*pi * gammaHz*Ga*ta;                              % [rad/cm]
dk = abs(k(2) - k(1));                                   % [rad/cm]

A    = exp(1i*(ones(M,1)*phi_e    + transpose(k)*x));     % [rad]    MxN (M lines, N columns)
A_HR = exp(1i*(ones(M,1)*phi_e_HR + transpose(k)*x_HR));  % [rad]    Mx(N*HRF) (M lines, N*HRF columns)

% ----------------------------------------
%  Perform the high-resolution evaluation
% ----------------------------------------
SR_modes = {'SR_iterative',...
	        'SR_inv'      ,...   % Manually calculate Matlab's pinv
			'SR_pinv'};
SR_mode = SR_modes{1};
switch SR_mode
	case 'SR_iterative'
		% (1)  Calculate an initial guess
		% - This might be done by FT-ing S, filtering all spectral components that are above
		%   the current resolution (meanining higher than 1/initial_dx) and IFT-ing the result.
% 		f0 = transpose(interp1(1:M,abs(S),linspace(1,M,N)));              % N
		% - Alternatively, just use zeros
		f0 = zeros(N,1);                                                  % N
		
		if (DEBUG_FLAG >= 3)
		figure; plot(linspace(x(1),x(end),M), abs(S)  / max(abs(S)) ,'k.-'); hold on;
		        plot(x                      , abs(f0) / max(abs(f0)),'r.-'); legend({'Original S','S0'});
		end;

		% (2)  Calculate weights around the stationary point, at each time-point
		for idx = 1:M
			v = A_HR(idx,:);                                              % N*HRF
			m = reshape(v,SR_HRF,N);                                         % NxHRF
			w_f = mean(m);                                                % N

			if (SRExp)
				[max_val,max_idx] = max(w_f);
				w_f = window_1D_vec_around_center_col(w_f,max_idx,SRExp,SRESP,0);
			end;

			W_f(idx,1:N) = w_f;                                           % MxN
		end;
		W_f = abs(W_f) / max(max(abs(W_f)));
		W_f = flipdim(W_f,2);
		
		[W_f_PSF, W_f_gauss] = calculate_spatio_temporal_weights(A_HR,N,M,SR_HRF,Ta,0,0,SRExp,SRESP,0,initial_dx/required_dx,dk,initial_dx,required_dx);
% 		W_f = W_f_PSF;
		
		if (DEBUG_FLAG >= 2)
% 		figure;imagesc(abs(W_f));
		A_tmp = ones(M,1)*phi_e_HR + transpose(k)*x_HR;
		A_tmp = flipdim(A_tmp,2);
		figure;
		subplot(1,2,1); imagesc(A_tmp);  set(gca,'YDir','normal'); title('A [rad]');
		subplot(1,2,2); imagesc(W_f);    set(gca,'YDir','normal'); title('W f [normalized]');
		end;
		A = A .* W_f;

		% (3)  Initialize variables
		y  = S;                                                           % M
		r0 = y - A*f0;                                                    % M = M - MxN * N
		W  = ones(M,1);  % Give the same weight to all data points        % M
		p0 = ctranspose(A) * (W.*r0);                                     % N = NxM * (M.*M)
		z0 = p0;                                                          % N
		n_iter = 1;
		
		% (4)  Iteratively solve for f
		f = increase_1Dvector_res(f0,S,A,M,x,n_iter);

% 		for idx = 1:n_iter
% 			v0 = A*p0;                                                    % M = MxN * N
% 			a0 = (ctranspose(z0)*z0) / (ctranspose(v0)*(W.*v0));          % 1 = 1/1 = (N' * N) / (M' * (M.*M))
% 			f  = f0 + a0*p0;                                              % N
% 			r  = r0 - a0*v0;                                              % M
% 			z  = ctranspose(A) * (W.*r);                                  % N = NxM * (M.*M)
% 			b0 = (ctranspose(z)*z) / (ctranspose(z0)*z0);                 % 1 = 1/1 = (N' * N) / (N' * N)
% 			p  = z + b0*p0;                                               % N
% 			
% 			% Set the next iteration values
% 			f0 = f;
% 			r0 = r;
% 			p0 = p;
% 			z0 = z;
% 		
% % 			figure; plot(linspace(x(1),x(end),M), abs(S) / max(abs(S)),'k.-'); hold on;
% % 			        plot(x                      , abs(f) / max(abs(f)),'r.-'); legend({'Original S','SR S'});
% % 			title(sprintf('SR profile (iter = %1.0f)',idx));
% 		end;
		new_S = f;
		
	case 'SR_inv'
		A_T = transpose(A);                                      % NxM
		A2 = A_T * A;                                            % NxN
		A2_inv = inv(A2);                                        % NxN
		new_S = A2_inv * A_T * S;                                % N    [NxN * (NXM * M) = NxN * (N) = N]
	
		disp(sprintf('Det = %3.3f\n',det(A2)));
		figure; plot(abs(eig(A2)),'.-'); title('A2 Eigen values');
		figure; imagesc(abs(A2_inv*A2)); title('(A^TA)^{-1}(A^TA) Matrix');

	case 'SR_pinv'
		A_pinv = pinv(A);                % NxM
		new_S = A_pinv * S;              % N    [NxM * M = N]
		
		disp(sprintf('Matrix rank = %3.3f\n',rank(A)));
		figure; imagesc(abs(A_pinv*A)); title('A^{-1}A Matrix');
		figure; imagesc(abs(A*A_pinv)); title('AA^{-1} Matrix');
end;

% Interpolate the new profile to the original (signal's) length
% new_S = interp1(1:N,new_S,linspace(1,N,M));

if (DEBUG_FLAG >= 2)
figure;
plot(abs(S)/max(abs(S)),'k.-'); hold on; plot(abs(new_S)/max(abs(new_S)),'r.-');
legend({'Initial S','New S'});
end;

return;
