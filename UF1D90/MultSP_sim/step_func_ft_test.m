clcl;
v = [zeros(1,64) ones(1,64)];
v = v .* exp(i*100*pi*v/length(v));
% xv = 1:length(v);
% pv = polyfit(xv,v,53);
% v = polyval(pv,xv);

for idx = 1:1
% 	v = transpose(smooth(v,5));
end;

ftv = fftshift(fft(fftshift(v)));
iftv = ifftshift(ifft(ifftshift(ftv)));

ftvf = window_1Dvec(ftv,20,5,1,'none');
iftvf = ifftshift(ifft(ifftshift(ftvf)));

figure;
subplot(1,3,1); plot(abs(v)    ,'.-'); title('orig v');
subplot(1,3,2); plot(abs(iftv) ,'.-'); title('iftv  ');
subplot(1,3,3); plot(abs(iftvf),'.-'); title('iftvf ');
