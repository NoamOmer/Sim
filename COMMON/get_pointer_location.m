function [xpos,ypos]=show_pointer_location;
x=get(0,'pointerlocation')+1;
a=get(gcf,'units');
b=get(gca,'units');
set(gcf,'units','pixels');
set(gca,'units','pixels');
p=get(gcf,'pos');
p1=get(gca,'pos');
xl=xlim;
yl=ylim;
xpos(1)=(x(1)-p(1)-p1(1))/p1(3)*diff(xl)+xl(1)+diff(xl)/2/p1(3);
ypos(1)=(x(2)-p(2)-p1(2))/p1(4)*diff(yl)+yl(1)+diff(yl)/2/p1(4);
xpos(2)=x(1)-p(1)+20;
ypos(2)=x(2)-p(2);
set(gcf,'units',a);
set(gca,'units',b);
txt=findobj(gcf,'tag','pointLoc');
if isempty(txt)
    txt=;
end
set(txt,'position',[x(2),y(2),80,20],'string',[num2str(x(1),4),' , ',num2str(y(1),4)]);
