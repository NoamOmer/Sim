% Window a 2D Matrix

% Inputs:
%  in_mat     : matrix : 2D matrix
%  center_row : scalar : Row   , around which the window should be centerd
%  center_col : scalar : Column, around which the window should be centerd
%  winF       : scalar : Window factor (for certain types of windows)
%  winType    : scalar : Window type
%                        (see MATLAB's help on 'window' function for a list of available windows)

function out_mat = window_2D_mat(in_mat,center_row,center_col,winF,winType)

[n_row,n_col] = size(in_mat);
% row_win = gausswin(2*n_row,gauss_factor*2);  % Use a more generalized window function
% col_win = gausswin(2*n_col,gauss_factor*2);
switch winType
case {'gausswin','kaiser','tukeywin'}
	row_win = window(winType,2*n_row,winF*2);
	col_win = window(winType,2*n_col,winF*2);
otherwise
	row_win = window(winType,2*n_row);
	col_win = window(winType,2*n_col);
end

win_mat = row_win * transpose(col_win);
win_mat = win_mat / max(max(win_mat));

 low_row = ((2*n_row)/2) - center_row;
high_row = low_row + n_row;
 low_col = ((2*n_col)/2) - center_col;
high_col = low_col + n_col;

% figure; imagesc(win_mat); set(gca,'YDir','normal');
% axis image;
win_mat = win_mat(low_row+1:high_row,low_col+1:high_col);
% figure; imagesc(win_mat); set(gca,'YDir','normal');
% axis image;
out_mat = in_mat .* win_mat;

return;

