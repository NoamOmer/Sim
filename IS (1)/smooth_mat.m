
function out_mat = smooth_mat(in_mat,SmoothFx,SmoothFy)

CONV_MAT = ones(SmoothFx,SmoothFy) / sum(sum(ones(SmoothFx,SmoothFy)));
out_mat = conv2(in_mat,CONV_MAT,'same');

return;

