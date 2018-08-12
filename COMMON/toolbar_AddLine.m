function toolbar_AddLine(a)
if nargin<1
    a=get(gcbo,'tag');
end
L=get(gcbo,'userdata');
l=[];
for k=1:length(L)
    if ishandle(L(k))
        l(end+1)=L(k);
    end
end
switch a
    case 'gcaHoriz'
        [x,y]=get_points(gcf,1);
        ax=get(gcf,'currentaxes');
        yl=get(ax,'ylim');
        l=line([x x],yl,'linestyle',':','linewidth',1.5,'color','r',...
            'buttondownfcn','linemotion(''moveX'')');
        set(l,'userdata',l);
    case 'Horiz'
        [x,y]=get_points(gcf,1);
        ax=findobj(gcf,'type','axes','visible','on');
        for k=1:length(ax)
            yl=get(ax(k),'ylim');
            l(k)=line([x x],yl,'linestyle',':','linewidth',1.5,'color','r','parent',ax(k));
        end
        S='linemotion(''moveX'',''toolbar_AddLine(''''AlignLines'''')'')';
        set(l,'userdata',l,'buttondownfcn',S);
    case 'AlignLines'
        x=get(gco,'xdata');
        set(l,'xdata',x)
        return
    case 'Vert'
        [x,y]=get_points(gcf,1);
        ax=get(gcf,'currentaxes');
        xl=get(ax,'xlim');
        l=line(xl,[y,y],'linestyle',':','linewidth',1.5,'color','r',...
            'buttondownfcn','linemotion(''moveY'')');
        set(l,'userdata',l);
    case 'Deleteline'
        set(0,'showhiddenhandles','on')
        delete(l)
        set(0,'showhiddenhandles','off')
        return
    case 'lineColor'
        set(0,'showhiddenhandles','on')
        c=uisetcolor(l(1));
        set(l,'color',c)
        set(0,'showhiddenhandles','off')
        return
    case 'lineStyle'
        set(0,'showhiddenhandles','on')
        q=get(l(1),'linestyle');
        s=set(l(1),'linestyle');
        k=find(strcmp(s,q)==1);
        [k,OK]=listdialog('liststring',s,'InitialValue',k,'SelectionMode','single','PromptString',...
            'Select line style','listsize',[100,100]);
        if OK==1
            set(l,'linestyle',s{k});
        end
        set(0,'showhiddenhandles','off')
        return
    otherwise
        errordlg('Unknown option');
        return
end


if length(l)==1
    for k=1:length(l)
        ucm=uicontextmenu;
        set(l(k),'handlevisibility','off','uicontextmenu',ucm)
        uimenu(ucm,'tag','Deleteline','userdata',l(k),...
            'label','Delete line','callback','toolbar_AddLine');
        uimenu(ucm,'tag','lineColor','userdata',l(k),...
            'label','Line color','callback','toolbar_AddLine');
        uimenu(ucm,'tag','lineStyle','userdata',l(k),...
            'label','Line style','callback','toolbar_AddLine');
    end
else
    for k=1:length(l)
        ucm=uicontextmenu;
        set(l(k),'handlevisibility','off','uicontextmenu',ucm)
        
        uimenu(ucm,'tag','Deleteline','userdata',l,...
            'label','Delete all lines','callback','toolbar_AddLine');
        
        uimenu(ucm,'tag','lineColor','userdata',l(k),...
            'label','Line color','callback','toolbar_AddLine');
        
        uimenu(ucm,'tag','lineStyle','userdata',l(k),...
            'label','Line style','callback','toolbar_AddLine');
        
%         uimenu(ucm,'tag','Deleteline','userdata',l,...
%             'label','Delete all lines','callback','toolbar_AddLine');
        
        uimenu(ucm,'tag','lineColor','userdata',l,...
            'label','Line color (all)','callback','toolbar_AddLine');
        
        uimenu(ucm,'tag','lineStyle','userdata',l,...
            'label','Line style (all)','callback','toolbar_AddLine');
    end    
end