function out_std = std_around_extremum(v);

N = length(v);
x = (1:N);
y = abs(v);

[yM,xM] = max(v);
y = y.*(y>(0.05*yM));   % Thresh-hold

tmp1 = (x.*y - xM*yM)./yM;
tmp2 = (1/(N-1)) * sum(tmp1.^2);
out_std = abs(xM - sqrt(tmp2));
