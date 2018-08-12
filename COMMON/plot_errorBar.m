function varargout=plot_errorBar(X,Y,Eup,varargin);
% Plots Y(X) with error bar
% Eup is the error as a function of X.
% In case of different positive and negative errors, apply the negative error vector as the 
% first value of varargin.
% The first output is a handle to the line Y(X) and the second is a handle
% to the line Error(X)
fig=gcf;
Edown=Eup;% setting the 
if length(varargin)>0
    if size(varargin{1})==size(Eup)
        Edown=varargin{1};
        varargin(1)=[];
    end
end
dx=mean(diff(X))/5;
X=X(:)';Y=Y(:)';Eup=Eup(:)';Edown=Edown(:)';
x=zeros(7,length(X));
x(1,:)=X-dx;   x(2,:)=X+dx;   x(3,:)=X;
x(4,:)=X;        x(5,:)=X-dx;    x(6,:)=X+dx;
x(7,:)=NaN*ones(1,length(X));
y=x;
for k=1:3,
    y(k,:)=Y-Edown;   
    y(k+3,:)=Y+Eup;      
end
y(7,:)=NaN*ones(1,length(Y));
l(1)=line(X,Y,'marker','O');
l(2)=line(x(:),y(:));
for k=1:2:length(varargin)
    set(l,varargin{k},varargin{k+1})
end
if nargout==1
    varargout{1}=l;
elseif nargout==2
    varargout{1}=l(1);
    varargout{2}=l(2);
end