function a=changePatchpattern(p,pattern,col,mark)
if nargin==0,return,end
d=30;
if findstr(lower(pattern),'full'),    d=50;end
ax=get(p,'parent');
xl=get(ax,'xlim');yl=get(ax,'ylim');
A=diff(yl)/diff(xl);
pos=get(ax,'pos');    
dy=diff(ylim)/d;
dx=diff(xlim)/d*pos(4)/pos(3);
x=get(p,'xdata');
y=get(p,'ydata');
switch pattern
    case {'FullNegLines','NegLines'}
        a=drawpattern(x,y,1,dy,dx,A);
    case {'FullPosLines','PosLines'}
        a=drawpattern(x,y,0,dy,dx,A);
    case {'FullCrossLines','CrossLines'}
        a=drawpattern(x,y,0,dy,dx,A);
        b=drawpattern(x,y,1,dy,dx,A);
        a=[a(:);b(:)];
    case {'FullSquareLines','SquareLines'}
        a=drawpattern(x,y,2,dy,dx,1);
    case {'FullDots','Dots'}
        [X,Y]=patchpattern(x,y,3,dy,dx,1);
        clear patchpattern
        a=line(X,Y,'tag','pattern','linestyle','none','marker',mark);        
end
set(a,'color',col,'hittest','off');

function a=drawpattern(x,y,pattern,dy,dx,A);
[X,Y]=patchpattern(x,y,pattern,dy,dx,A);
clear patchpattern
x1=vec2mat(X,2,0);
x1(3,:)=NaN;
y1=vec2mat(Y,2,0);
y1(3,:)=NaN;
a=line(x1(:),y1(:),'tag','pattern');
