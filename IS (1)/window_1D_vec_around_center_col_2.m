% Window a 1D vector

% Inputs:
%  in_vec     : vector : 1D vector
%  center_col : scalar : Column, around which the window should be centerd
%  winF       : scalar : Window factor (for certain types of windows)
%  winType    : scalar : Window type
%                        (see MATLAB's help on 'window' function for a list of available windows)
function out_vec = window_1D_vec_around_center_col_2(in_vec,center_col,winF,winType)

n_col = length(in_vec);
switch winType
case {'gausswin','kaiser','tukeywin'}
	col_win = transpose(window(winType,2*n_col,winF*2));
otherwise
	col_win = transpose(window(winType,2*n_col));
end

col_win = col_win / max(max(col_win));

 low_col = ((2*n_col)/2) - center_col;
high_col = low_col + n_col;

% figure; imagesc(col_win); set(gca,'YDir','normal'); axis image; hold on;
col_win = col_win(low_col+1:high_col);
% figure; imagesc(col_win); set(gca,'YDir','normal'); axis image;

out_vec = in_vec .* col_win;

return;

