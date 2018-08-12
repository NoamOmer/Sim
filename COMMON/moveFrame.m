% function moveFrame(h,dx,dy)
% f=get(h,'parent');
% a=get(f,'children');
% pPrev=get(h,'pos');
% for k=1:length(a)
%     P=get(a(k),'pos');
%     if P(1)>=pPrev(1) & P(1)+P(3)<=pPrev(1)+pPrev(3) & P(2)>=pPrev(2) & P(2)+P(4)<=pPrev(2)+pPrev(4)
%         P(1)=P(1)+dx;P(2)=P(2)+dy;
%         set(a(k),'pos',P)
%     end
% end
% 
function moveFrame(opt)
h=gco;
switch opt
    case 'initiate'
        f=get(h,'parent');
        a=get(f,'children');
        ph=get(h,'pos');
        l=1;
        for k=1:length(a)
            P=get(a(k),'pos');
            if P(1)>=ph(1) & P(1)+P(3)<=ph(1)+ph(3) & P(2)>=ph(2) & P(2)+P(4)<=ph(2)+ph(4)
                u(l)=a(k);
                pos(l,:)=P;
                l=l+1;
            end
        end
        X=get(gcf,'currentpoint');
        S.x=X(1);
        S.y=X(2);
        S.hand=u;
        S.position=pos;
        set(h,'userdata',S)
        s=get(gcf,'windowbuttonupfcn');        
        set(gcf,'doublebuffer','on','windowbuttonupfcn',['moveframe(''released'');%',s]);
        s=get(gcf,'windowbuttonmotionfcn');        
        set(gcf,'pointer','cross','windowbuttonmotionfcn',['moveframe(''moveframe'');%',s]);
    case 'released'
        s=get(gcf,'windowbuttonupfcn');
        l=findstr(s,'%');s(1:l(1))=[];
        set(gcf,'pointer','arrow','windowbuttonupfcn',s)   
        eval(s);
        s=get(gcf,'windowbuttonmotionfcn');
        l=findstr(s,'%');s(1:l(1))=[];
        set(gcf,'pointer','arrow','windowbuttonmotionfcn',s)
    case 'moveframe'
        S=get(h,'userdata');
        X=get(gcf,'currentpoint');
        x=S.x;
        y=S.y;
        u=S.hand;
        p=S.position;
        dx=X(1)-x;
        dy=X(2)-y;
        p(:,1)=p(:,1)+dx;        
        p(:,2)=p(:,2)+dy;
        for k=1:length(u)
            set(u(k),'pos',p(k,:))
        end
        setfield(S,'x',X(1));
        setfield(S,'y',X(2));
        set(h,'userdata',S)
end