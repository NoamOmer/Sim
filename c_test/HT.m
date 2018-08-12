function y=HT(x)
%Hilbert transform

M=length(x);
k=0:M-1;
s=zeros(1,M);

%s=s+((k>0) & (k<M/2))-((k>M/2) & (k<M));

s=s+(k==0)+(k==M/2)+2*((k>0) & (k<M/2));

y=imag(fft(s.*ifft(x)));

