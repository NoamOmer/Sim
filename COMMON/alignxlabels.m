function alignxlabels(fig,pos)
if nargin<1
    fig=gcf;
end
ax=findobj(fig,'type','axes')
for k=1:length(ax)
    q(k)=get(ax(k),'xlabel');
    set(q(k),'units','norm');
    p(k,1:3)=get(q(k),'position');
end

if nargin>1
    set(q,'position',pos)
else
    x=0.5;
    y=min(p(:,1));
    z=0;
    set(q,'position',[x,y,z])
end
