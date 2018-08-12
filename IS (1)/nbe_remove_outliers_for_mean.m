% Remove points that reside outside mean +- N x STD of original data
% Remove illegal values
function out_vec = nbe_remove_outliers_for_mean(vec,ignored_val,min_val,max_val,nSTD)
if (isempty(vec))
	out_vec=vec;
	return;
end;

out_vec = vec(vec~=ignored_val);
out_vec = out_vec(out_vec>min_val);
out_vec = out_vec(out_vec<max_val);

vec_mean = mean(out_vec);
vec_std  =  std(out_vec);

if (nSTD>0)
out_vec(out_vec < max([0,(vec_mean - nSTD*vec_std)])) = ignored_val;
out_vec(out_vec >        (vec_mean + nSTD*vec_std) ) = ignored_val;
out_vec = out_vec(out_vec~=ignored_val);
end;

return;
