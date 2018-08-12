
function [padded_mat n_rows n_cols] = zero_pad_1D(mat)

[n_rows ,n_cols] = size(mat);
two_powers = 2.^(1:35);

% Transpose vector, if it is not in a column format
if (n_cols == 1)
    mat = transpose(mat);
    [n_rows ,n_cols] = size(mat);
end;

% if (n_rows == 1)
% 	% Vector
	loc = find(two_powers > n_cols);
	pad_sz = two_powers(loc(1)) - n_cols;

	pad_left_sz   = round(pad_sz/2);
	pad_right_sz  = pad_sz - pad_left_sz;
	pad_left_vec  = zeros(n_rows,pad_left_sz);
	pad_right_vec = zeros(n_rows,pad_right_sz);

	padded_mat = [pad_left_vec mat pad_right_vec];
	[n_rows ,n_cols] = size(padded_mat);
% else
% 	error('>1D not supported ... use zero_pad_2D function instead');
% end;

return;
