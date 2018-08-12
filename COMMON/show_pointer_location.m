function show_pointer_location(opt);
txt=findobj(gcf,'tag','pointLoc');
if isempty(txt)
    txt=add_button('text',[0 0 .01,.01],'',[],'pointLoc');
    set(txt,'units','pixel','fontname','Times New Roman','fontsize',8,'horizontalalignment','left',...
        'foregroundcolor','r')
end

switch opt
case 1 % Show tracks
    units0=get(0,'units');
    unitsGcf=get(gcf,'units');
    
    set([0,gcf],'units','pixels');

    x=get(0,'pointerlocation')+1;
    p=get(gcf,'pos');
    ax=findobj(gcf,'type','axes','visible','on');
    flag=0;
    for k=1:length(ax)
        unitsGca=get(ax(k),'units');
        set(ax(k),'units','pixels');
        p1=get(ax(k),'pos');
        xl=get(ax(k),'xlim');
        yl=get(ax(k),'ylim');
        xpos=(x(1)-p(1)-p1(1))/p1(3)*diff(xl)+xl(1)+diff(xl)/2/p1(3);
        ypos=(x(2)-p(2)-p1(2))/p1(4)*diff(yl)+yl(1)+diff(yl)/2/p1(4);
        set(ax(k),'units',unitsGca);
        if prod(xpos-xl)<0 &prod(ypos-yl)<0
            flag=1;
            break
        end
    end
    set(0,'units',units0)
    set(gcf,'units',unitsGcf);
    if flag==0,delete(txt),return;end
    
    set(txt,'string',{num2str(xpos,4);num2str(ypos,4)},'visible','on');

    y=x(2)-p(2)-30;
    x=x(1)-p(1)+20;
    
    pTxt=get(txt,'extent');
    p1=pTxt+[x,y,0,0];
    if p1(1)+p1(3)>p(3),
       p1(1)=p1(1)-p1(3)-30;
   end
   if p1(2)<0,
       p1(2)=y+20;
   end
   set(txt,'position',p1)
case 0 % Do not track
    delete(txt)
    s=get(gcf,'windowbuttonmotion');
    k=findstr(lower(s),'show_pointer_location(0);');
    s(k+(0:24))=[];
    set(gcf,'windowbuttonmotion',s)
end
