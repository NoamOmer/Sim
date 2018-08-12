function [row_idx,col_idx] = max_mat(mat)
[col_max_vals,col_max_idcs] = max(mat);
[row_max_val ,row_max_idx ] = max(col_max_vals);

col_idx = row_max_idx;
row_idx = col_max_idcs(row_max_idx);

return;
