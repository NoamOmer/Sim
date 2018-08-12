function Linemotion(opt,action,h)
set(0,'showhiddenhandles','on')
if nargin==0,    opt='moveXY'; end
if nargin<2,      action =''; end
if nargin<3,      h=gco;end

switch opt
    case {'moveXY','moveX','moveY'}
        s=get(gcf,'windowbuttonupfcn');        
        set(gcf,'doublebuffer','on','windowbuttonupfcn',...
            ['linemotion(''released'');',action,';%',s]);
        x=get(gca,'currentpoint');        
        preUD=get(h,'userdata');
        if strcmp(opt,'moveY'),D.x=[];else,D.x=x(1,1);end    
        if strcmp(opt,'moveX'),D.y=[];else,D.y=x(1,2);end
        D.preUD=preUD;
        set(h,'userdata',D);
        s=get(gcf,'windowbuttonmotionfcn');    
        S=['linemotion(''moveline'');;;',s];
        if strcmp(lower(get(h,'type')),'patch')
            if ~isempty(preUD)
                if ishandle(preUD(1)),
                    S=['linemotion(''movePatchStruct'');;;',s];
                end
            end
        end
        if  isfield(get(h),'Position')
            S=['linemotion(''moveShape'');;;',s];
        end
        set(gcf,'pointer','cross','windowbuttonmotionfcn',S);
    case 'released'
        a=get(gco,'userdata');
        set(gco,'userdata',a.preUD)
        s=get(gcf,'windowbuttonupfcn');
        S='linemotion(''released'');';
        l1=findstr(s,S);
        s(l1-1+(1:length(S)))=[];

        l=findstr(s,'%');
        action=s(1:l(1));
        s(1:l(1))=[];
        set(gcf,'pointer','arrow','windowbuttonupfcn',s)   
        eval(action);
        s=get(gcf,'windowbuttonmotionfcn');
        l=findstr(s,';;;');s(1:(l(1)+2))=[];
        set(gcf,'pointer','arrow','windowbuttonmotionfcn',s)
        set(0,'showhiddenhandles','off')
    case 'moveline'
        a=get(h(1),'userdata');
        x=get(gca,'currentpoint');x=x(1,1:2);
        X=get(h,'xdata');
        Y=get(h,'ydata');
        xx=a.x;if ~isempty(xx), set(h,'xdata',X-xx+x(1));       a.x=x(1,1);end
        yy=a.y;if ~isempty(yy), set(h,'ydata',Y-yy+x(2));       a.y=x(1,2);end        
        set(h,'userdata',a)
    case 'movePatchStruct'
        a=get(h,'userdata');
        x=get(gca,'currentpoint');x=x(1,1:2);
        ha=a.preUD;
        for k=1:length(ha)
            X=get(ha(k),'xdata');
            Y=get(ha(k),'ydata');
            xx=a.x;if ~isempty(xx), set(ha(k),'xdata',X-xx+x(1));end
            yy=a.y;if ~isempty(yy), set(ha(k),'ydata',Y-yy+x(2));end
        end
        X=get(h,'xdata');
        Y=get(h,'ydata');
        xx=a.x;if ~isempty(xx), set(h,'xdata',X-xx+x(1));a.x=x(1,1);end
        yy=a.y;if ~isempty(yy), set(h,'ydata',Y-yy+x(2));a.y=x(1,2);end
        set(h,'userdata',a)
    case 'moveShape'
        a=get(gco,'userdata');
        x=get(gca,'currentpoint');x=x(1,1:2);
        X=get(gco,'position');
        xx=a.x;if ~isempty(xx), X(1)=X(1)-xx+x(1);       a.x=x(1,1);end
        yy=a.y;if ~isempty(yy), X(2)=X(2)-yy+x(2);       a.y=x(1,2);end
        set(gco,'userdata',a,'position',X)
    otherwise 
        errordlg(['Not a recognized option (opt=',opt,')'])
end