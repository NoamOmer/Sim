function l=plot_colored_line(x,y,c);
if nargin==1
    y=x;
    x=1:length(y);
    c=colormap;
end
if nargin==2
    c=colormap;
end
[minY,maxY]=min_max(y);
rangeY=maxY-minY;
K=round((y-minY)/rangeY*(length(c)-1))+1;
l=[];
for k=2:length(y)
    l(k-1)=line([x(k-1),x(k)],[y(k-1),y(k)],'color',c(K(k),:));
end