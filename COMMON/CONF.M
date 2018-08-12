function conf(x,y,z,levels),
if nargin<4
   levels=10;
end

contourf(x,y,z,levels)
c=findobj(gca,'type','patch');
cc=get(c,'cdata');
ca=caxis;
cl=cat(1,cc{:});
l=length(colormap);
I=round((cl-ca(1))/(ca(2)-ca(1))*(l-1))+1;
map=colormap;
for q=1:length(c)
   if ~isnan(cl(q))
      set(c(q),'edgecolor',map(I(q),:))
   else
      set(c(q),'edgecolor',[ 1 1 1]);
   end
end
