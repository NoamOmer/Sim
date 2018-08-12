% M   number of time points
% N   Number of pixels (required from the SR) on the PE axis
function [W_f, W_f_gauss] = calculate_spatio_temporal_weights(A_HR,N,M,HRF,Ta,T2PE_factor,plot_flag,SRExp,SRESP,SR_ones_win_flag,SR_factor,dk_rad,dz_LR,dz_HR)
set_globals;
dk = dk_rad / (2*pi);
t  = linspace(0,Ta,M);
if (T2PE_factor)
	t2_relax = exp(-t/T2PE_factor);
else
	t2_relax = ones(1,M);
end;

for idx = 1:M
	v = A_HR(idx,:);                                                  % N*HRF
	m = reshape(v,HRF,N);                                             % NxHRF
	
	% Calculate the intra pixel dephasing
	q=sqrt(N);
	if (1)
	exp_arg = ( idx*dk*(dz_LR^2)  +  (-N:-1)*dz_HR ) / (dz_LR*2*pi*q);
	w_f_gauss = exp(-0.5 * (exp_arg.^2) );                            % [nbe oct 21 2011] fit to theory
	else
	w_f_gauss = exp(-0.5 * ((((1:M) - idx + 0.5) / SRExp).^2) );      % [nbe oct 13 2011] unchanged
	w_f_gauss = w_f_gauss(round(linspace(1,M,N)));                    % [nbe oct 13 2011] added
	w_f_gauss = w_f_gauss(end:-1:1);
	end;
	
	w_f = abs(mean(m));                                               % N

% 	w_f = transpose(smooth(w_f,20)); % smooth: causes degradation of the SR process
% 	w_f = w_f.^1.9;                  % window: another indirect method for windowing the weighting vector

	[max_val,max_idx] = max(w_f);    % find the current maximal weight

	if (SR_ones_win_flag)
		w_f = ones(1,N);             % give identical weight to all points - this must be accompanied with SRExp~=0
	end;
% 	if (SRExp == 0 && SR_ones_win_flag), uiwait(msgbox('Warning: Using equal weight matrix w/o windowing')); error('Exiting'); end;

	% Limit the weighting response function around the SP region
	if (SRExp)
% 	w_f       = w_f .* w_f_gauss; % [nbe 2012_03_30 - replace windowing function with a simple (yet effective) gaussian]
	w_f       = window_1D_vec_around_center_col(w_f      ,max_idx,SRExp,SRESP,0);     % window around the incremented index
	w_f_gauss = window_1D_vec_around_center_col(w_f_gauss,max_idx,SRExp,SRESP,0);     % window around the incremented index
	end;
	W_f(idx,1:N) = w_f;                                               % MxN
	W_f_gauss(idx,1:N) = w_f_gauss;
	
	% T2 correct the amplitude
	W_f(idx,1:N)       = W_f(idx,1:N)       * t2_relax(M-idx+1);
	W_f_gauss(idx,1:N) = W_f_gauss(idx,1:N) * t2_relax(M-idx+1);
end;

% Normalize the weighting function to [0..1]
W_f = abs(W_f) / max(max(abs(W_f)));
W_f_gauss = abs(W_f_gauss) / max(max(abs(W_f_gauss)));

if (plot_flag)  figure;imagesc(abs(W_f));
                figure;imagesc(abs(W_f_gauss));  end;

return;

