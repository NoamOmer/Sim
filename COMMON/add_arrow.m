function [ArrowLine,ArrowHead]=add_arrow(varargin)
% This function adds an arrow into the current axis
% (you may call the function without arguments)
% the first 2 input arguments are X and Y, which are the locations of the head and tail of the arrow
% additional input should be pairs of property/value for the arrow.
% for example: [a,b]=add_arrow([5 5],[3,4],'color','g','facecolor','r');
% the output arguments are handles to the line and head of the arrow
if nargin==0
    [X,Y]=get_points(gca,2);
else 
    X=varargin{1};Y=varargin{2};
    varargin(1:2)=[];
end

ArrowLine=line(X,Y,'linewidth',2,'linestyle','-','color','K','tag','arrowline');

p=get(gca,'pos');
p1=get(gcf,'pos');
dx=p(3)*p1(3);
dy=p(4)*p1(4);

rx=diff(xlim)/40/dx/2;
ry=diff(ylim)/30/dy/2;
alpha=atan2(diff(Y)/ry,diff(X)/rx);

c=abs(cos(alpha));
s=abs(sin(alpha));
sx=sign(diff(X));
sy=sign(diff(Y));
if sx==0,sx=1;end
if sy==0,sy=1;end

x=X(2)+rx*sx*[-s/2-c,0,s/2-c,-s/2-c];
y=Y(2)+ry*sy*[c/2-s,0,-c/2-s,c/2-s];

ArrowHead=patch(x,y,'K');
set(ArrowHead,'edgecolor','K','tag','arrowhead')
head_prop=lower(fieldnames(get(ArrowHead)));
line_prop=lower(fieldnames(get(ArrowLine)));
for k=1:2:length(varargin)
    if sum(strcmp(line_prop,varargin{k}))>0
        set(ArrowLine,varargin{k},varargin{k+1})
    end
    if sum(strcmp(head_prop,varargin{k}))>0
        set(ArrowHead,varargin{k},varargin{k+1})
    end
end