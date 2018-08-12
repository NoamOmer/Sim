% Window a 1D vector using an EXPONENTIAL window
% Assumptions:
%  1. vector is either positive or negative
%  2. vector is real

% Inputs:
%  v       : vector : 1D vector
%  exp_val : scalar : Windowing exponent value
%  ESP     : scalar : Edge slope parameter
function [v_out] = window_1D_vec(v_in,exp_val,ESP,DEBUG_FLAG,title_str)
v = v_in;
temp = v(round(length(v)/2));
if (imag(temp) == 0)
	sv = sign(v(round(length(v)/2)));
else
	sv = 1;
end;
v  = abs(v);
mv = max(v);
v  = v/mv;
tmp1 = length(v);
tmp2 = ((-round(tmp1/2)+1):1:(+round(tmp1/2)));
win = exp(-(exp_val*abs(tmp2./tmp1)).^ESP);
win = win(1:tmp1);

windowed_v = v .* win;
v_out = v_in .* (1e-10 + win/max(win));
windowed_v = windowed_v * mv;   v = v*mv; % return the original scaling
windowed_v = windowed_v * sv;   v = v*sv; % return the original sign
if (DEBUG_FLAG >= 1)
	figure; hold;
	plot(1:length(v), v         ,'b.-');
	plot(1:length(v), windowed_v,'g.-');
	title(title_str);
	legend({'pre-win','post-win'},'Location','Best');
end;

return;

