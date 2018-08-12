
function [padded_mat_3D n_rows n_cols n_planes nFourthDim] = zero_pad_2D(mat)

orig_mat_3D = mat;
[n_rows ,n_cols, n_planes, n_fourth_Dim] = size(orig_mat_3D);
two_powers = 2.^(1:35);

for forthDim_idx = 1:n_fourth_Dim
	for plane_idx = 1:n_planes
		clear mat;
		mat = orig_mat_3D(:,:,plane_idx);

		% 1. pad columns
		loc = find(two_powers > n_cols);
		pad_sz        = two_powers(loc(1)) - n_cols;
		pad_left_sz   = round(pad_sz/2);
		pad_right_sz  = pad_sz - pad_left_sz;
		pad_left_mat  = zeros(n_rows,pad_left_sz);
		pad_right_mat = zeros(n_rows,pad_right_sz);

		padded_mat = [pad_left_mat mat pad_right_mat];

		% 2. pad rows
		loc = find(two_powers > n_rows);
		pad_sz       = two_powers(loc(1)) - n_rows;
		pad_up_sz    = round(pad_sz/2);
		pad_down_sz  = pad_sz - pad_up_sz;
		pad_up_mat   = zeros(pad_up_sz,size(padded_mat,2));
		pad_down_mat = zeros(pad_down_sz,size(padded_mat,2));

		padded_mat = [pad_up_mat; padded_mat; pad_down_mat];
		padded_mat_3D(:,:,plane_idx,forthDim_idx) = padded_mat;
	end;
end;

[n_rows, n_cols, n_planes, n_fourth_Dim] = size(padded_mat);

return;
