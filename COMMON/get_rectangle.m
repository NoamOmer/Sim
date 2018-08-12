function R=get_rectangle(opt)
if nargin==1
    feval(opt)
    return
end
a=findobj(gcf,'type','axes','visible','on');
if isempty(a),a=gca;end
set(0,'showhiddenhandles','on')
hitTest=findobj(gcf,'hittest','on');
set(hitTest,'hittest','off')
G.XlimMode=get(gca,'xlimmode');
G.YlimMode=get(gca,'ylimmode');
set(a,'hittest','on','ylimmode','manual','xlimmode','manual')
uistate=uisuspend(gcf);
set_pointer('getPoints');
S=get(gcf,'userdata');
A=get(gcf,'windowbuttonmotionfcn');
for k=1:length(a)
    U{k}=get(a(k),'buttondownfcn');
end
set(gcf,'userdata',[],'doublebuffer','on')
set(a,'buttondownfcn','get_rectangle(''buttonDown'');')
waitfor(gcf,'userdata','completed')
set(gcf,'windowbuttonmotionfcn',A)
set(gcf,'userdata',S)
for k=1:length(a)
    set(a(k),'buttondownfcn',U{k});
end
r=findobj(gcf,'tag','temp_rect');
x=get(r,'xdata');mx=min_max(x);
y=get(r,'ydata');my=min_max(y);
R=[mx(1),my(1),diff(mx),diff(my)];
delete(r)
uirestore(uistate);
set(hitTest,'hittest','on')
set(a,'ylimmode',G.YlimMode,'xlimmode',G.XlimMode)
set(0,'showhiddenhandles','off')
function buttonDown
X=get(gca,'currentpoint');X=[X(1,1),X(1,2)];
dx=diff(xlim)/100;dy=diff(ylim)/100;
x=X(1)+[0,0,dx,dx,0];
y=X(2)+[0,dy,dy,0,0];
r=line(x,y,'color','r','linestyle',':','tag','temp_rect');
set(gcf,'userdata',[r,X],'windowbuttonmotionfcn','get_rectangle(''buttonMotion'');',...
    'windowbuttonupfcn','get_rectangle(''buttonUp'');')

function buttonMotion
X=get(gca,'currentpoint');X=[X(1,1),X(1,2)];
U=get(gcf,'userdata');
r=U(1);
x=U(2);
y=U(3);
set(r,'xdata',[x,x,X(1),X(1),x],'ydata',[y,X(2),X(2),y,y]);

function buttonUp
set(gcf,'userdata','completed')
