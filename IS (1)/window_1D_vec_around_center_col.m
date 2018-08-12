% Window a 1D vector around a center column

% Inputs:
%  in_vec     : vector : 1D vector
%  center_col : scalar : Column, around which the window should be centerd
%  exp_val    : scalar : Windowing factor
%  ESP        : scalar : Windowing Edge Slope Parameter

function out_vec = window_1D_vec_around_center_col(in_vec,center_col,exp_val,ESP,plot_flag,ax)

v = in_vec;
temp = v(round(length(v)/2));
if (imag(temp) == 0)
	sv = sign(v(round(length(v)/2)));
else
	sv = 1;
end;
v  = abs(v);
mv = max(v);
v  = v/mv;
tmp1 = length(v)*2;  % Create a window with twice the length of the vector 
tmp2 = ((-round(tmp1/2)+1):1:(+round(tmp1/2)));
win = exp(-(exp_val*abs(tmp2./tmp1)).^ESP);
% win = win(1:tmp1);

n_col    = length(in_vec);
low_col  = ((2*n_col)/2) - center_col;
high_col = low_col + n_col;
col_win  = win(low_col+1:high_col);

out_vec = in_vec .* col_win;
fh = 0;
if (plot_flag)
	if (exist('ax','var') ~= 0) && (ax ~= 0),
% 		delete(get(ax,'Children'));
		axes(ax);
		hold off;
	else
		fh = figure;
	end;
% 	if (fh)
	plot(abs(in_vec ) / max(abs(in_vec )),'b.-'); hold on;
	plot(abs(out_vec) / max(abs(out_vec)),'g.-');
	plot(abs(col_win) / max(abs(col_win)),'r'  ); % legend({'orig vec','filtered vec','window'});
% 	end;
end;
return;

