function alignYlabels(fig,pos)
if nargin<1
    fig=gcf;
end
ax=findobj(fig,'type','axes')
for k=1:length(ax)
    q(k)=get(ax(k),'ylabel');
    set(q(k),'units','norm');
    p(k,1:3)=get(q(k),'position');
end

if nargin>1
    set(q,'position',pos)
else
    x=min(p(:,1));
    y=0.5;
    z=0;
    set(q,'position',[x,y,z])
end
