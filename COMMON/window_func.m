function w=window(varargin)
N=varargin{1};WIN=lower(varargin{2});
if length(varargin)==3
   beta=varargin{3};
end

switch WIN
case 'rectangle',w=ones(N,1);
case 'hanning',w=hanning(N);
case 'hamming',w=hamming(N);
case 'blackman',w=blackman(N);
case 'kaiser',w=kaiser(N,beta);
case 'boxcar',w=boxcar(N);
case 'square',w=boxcar(N);
case 'triangle',w=bartlett(N);
case 'bartlett',w=bartlett(N);
case 'rao',w=rao(N,beta);  
case 'buchman',
   m=0:(N/2-1);
   w=1/pi*abs(sin(2*pi*m/N))+(1-2*abs(m)/N).*cos(2*pi*m/N);
   w=[w(end:-1:1),w]';   % '
end
w=w/sum(w);

function opwind=rao(nfft,winsize)
winsize = winsize - rem(winsize,2) + 1;  % make it odd
if (winsize > 1)
   mwind   = fix (nfft/winsize);            % the scale parameter M
   lby2    = (winsize - 1)/2;
   theta  = -lby2:lby2;
   opwind = ones(winsize,1) * (theta .^2);       % w(m,n)=m^2
   opwind = opwind + opwind' + theta' * theta;   % m^2 + n^2 + mn
   opwind = 1 - (2*mwind/nfft)^2 * opwind ;      %
   hex    = ones(winsize,1) * theta;             % m
   hex    = abs(hex) + abs(hex') + abs(hex+hex');
   hex    = (hex < winsize);
   opwind = opwind .* hex;
   opwind = opwind * (4 * mwind^2) / (7 * pi^2) ;
end
opwind=opwind/sum(opwind(:));


