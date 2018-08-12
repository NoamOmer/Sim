function  [x1,x2,x3,x4]=kk(n,Y1,Y2,Y3,Y4)
% The function get up to 4 vectors and return the linear
% detrending of every n points.
m=nargin-1;
if m==1;Y2=0;Y3=0;Y4=0;end
if m==2;Y3=0;Y4=0;end
if m==3;Y4=0;end
if nargin~=(nargout+1)
   disp('please enter the number of points interval ');
   disp('or the apropriate number of input/output vectors');return
end
if max(size(n))>1
   a=[1;n];
else
   a=[1:n:length(Y1) length(Y1)];
end
if (min(size(Y1))>1 |min(size(Y2))>1 |min(size(Y3))>1 |min(size(Y4))>1)
   error('expecting vectors (not matrics) ');return
end

nargout=m;
x1=[];x2=[];x3=[];x4=[];
for mm=1:m
    C=['X=Y',num2str(mm),';'];
    eval(C);
    x=[];
    for k=1:(length(a)-1)   
       l=a(k):(a(k+1));
       M=detrend(X(l));
       x(l)=M;
    end
    C=['x',num2str(mm),'=x;'];
    eval(C);
 end