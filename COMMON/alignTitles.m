function alignTitles(fig,pos,varargin)
if nargin<1
    fig=gcf;
end
ax=findobj(fig,'type','axes');
for k=1:length(ax)
    q(k)=get(ax(k),'title');
    set(q(k),'units','normalized');
    p(k,1:3)=get(q(k),'position');
end

if nargin>1 & length(pos)==3
    set(q,'position',pos)
else
    x=0.5;
    y=1.01;
    z=0;
    set(q,'position',[x,y,z])
end
for k=1:2:length(varargin)
    set(q,varargin{k},varargin{k+1})
end