
% Assumption: all w are non-zero
function [wstd_val] = wstd(v,w,wmean_val)

N=length(v);
if ~exist('wmean_val','var')
	wmean_val = wmean('arithmetic',v,w,2);
end;
wstd_val = sqrt( sum(w.*((v-wmean_val).^2)) / ...
             ((N-1) * sum(w) / N) );
return;

