function out_std = std2_around_extremum(v);

N = length(v);
x = (1:N);
y = abs(v);

[yM,xM] = max(v);
y = y.*(y>(0.05*yM));   % Thresh-hold

tmp1 = sum(((x-xM).^2) .* y);
tmp2 = sum(y);
out_std = sqrt(tmp1/tmp2);
