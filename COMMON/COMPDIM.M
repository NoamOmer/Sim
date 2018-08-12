function [C,d,R,H,B,r]=compdim(x,tau,embed_dim,rmin,rmax,num_bin,skip_factor,flag);
% [C,d,R,h]=compdim(x,tau,embed_dim,rmin,rmax,num_bin,skip_factor,flag);
% Log of the Correlation sum
% d - embedding dimensions
% R - Radius
%
% x- The signal
% tau - tau
% embed_dim - maximal embedding dimension
% rmin - minimal distance
% rmax - maximal distance
% num_bin - number of bins
% skip factor- skip factor

if nargin<7
   errror('You must enter 7 input arguments');
end
if nargin<8
   flag='fast';
end

tic;
if ~strcmp(flag,'fast')
   [H,d,r]=gpdim(x,tau,embed_dim,rmin,rmax,num_bin,skip_factor);
else
   [H,d,r]=gpdim_simple(x,tau,embed_dim,rmin,rmax,num_bin,skip_factor);
end
[m,n]=size(H);
R=zeros(size(H));
b=hanning(21);b=b/sum(b);
h=H/sum(H(:,1));
h(1,:)=zeros(1,size(h,2));
h=filtfilt(b,1,h);
C=log10(cumsum(h));
for k=1:n
   R(:,k)=r-log10(d(k));
end
B=diff(C)/(r(2,1)-r(1,1));
r=R(2:size(R,1),:);
y=toc;
disp(['CPU time: ' num2str(y)])
