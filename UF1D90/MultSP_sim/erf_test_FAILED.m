dx = 0.1;
x  = -3:dx:3;
erfx = (2/sqrt(pi)) * cumsum(exp(-(x.^2)))*dx - 1;
% erfx = erf(x);
figure;
subplot(1,2,1); plot(erfx,'.-'); title('erf(x)');

z  = 1*x;
n = 1:10000;
erfz = 0;
for idx = 1:length(n)
	n_series = ((-1)^n(idx)) .* (z.^(2*n(idx)+1)) ./ (factorial(n(idx)) * (2*n(idx)+1));
	erfz = erfz + n_series;
end;
erfz = erfz * (2/sqrt(pi));
subplot(1,2,2); plot(erfz,'.-'); title('erf(z)');