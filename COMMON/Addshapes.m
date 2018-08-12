function l=Addshapes(varargin)
if nargin==0
    a=get(gcbo,'tag');
    col='r';
elseif isstr(varargin{1})
    a=varargin{1};
    col='r';
else
    a='Manual';
end
switch a
case 'rectangle'
    R=get_rectangle;
    x=R(1)+[0,0,R(3),R(3),0];
    y=R(2)+[0,R(4),R(4),0,0];
case 'triangle'
    R=get_rectangle;
    x=R(1)+[0,R(3)/2,R(3),0];
    y=R(2)+[0,R(4),0,0];
case 'ellipse'        
    R=get_rectangle;
    a=R(3)/2;b=R(4)/2;
    x=linspace(-a,a,200);
    y=b/a*sqrt(a*a-x.*x);
    x=[x,fliplr(x)];y=[y,-y];
    x(end+1)=x(1);y(end+1)=y(1);
    x=x+a+R(1);y=y+b+R(2);
case 'Manual'
    x=varargin{1};y=varargin{2};col=varargin{3};
case 'Deleteline'
    set(0,'showhiddenhandles','on')
    l=get_objects;
    delete(l)
    set(0,'showhiddenhandles','off')
    return
case 'Edit'
    set(0,'showhiddenhandles','on')
    l=get(gcbo,'userdata');
    patcheditor(l);
    return
case 'backward'
    set(0,'show','on')
    l=get_objects;
    ax=get(l(1),'parent');
    a=get(ax,'children');
    for m=1:length(l)
        k=find(a==l(m));
        a(k)=[];
        a(end+1)=l(m);
    end
    set(ax,'children',a);
    set(0,'show','off')
    return
case 'front'
    set(0,'show','on')
    l=get_objects;
    ax=get(l(1),'parent');
    a=get(ax,'children');
    for m=1:length(l)
        k=find(a==l(m));
        a(k)=[];
        a(end+1)=l(m);
        a=[l(m);a(:)];
    end
    set(ax,'children',a);
    set(0,'show','off')
    return    
otherwise
    errordlg('Unknown option');
    return
end
l=patch(x,y,col);
ucm=uicontextmenu;
set(l,'linestyle','-','linewidth',1.5,'edgecolor','r','buttondownfcn','linemotion(''moveXY'')',...
    'handlevisibility','off','uicontextmenu',ucm);

uimenu(ucm,'tag','Deleteline','userdata',l,...
    'label','Delete object','callback','Addshapes');
uimenu(ucm,'tag','Edit','userdata',l,...
    'label','Edit shape','callback','Addshapes');
uimenu(ucm,'tag','backward','userdata',l,...
    'label','Send to back','callback','Addshapes');
uimenu(ucm,'tag','front','userdata',l,...
    'label','Bring to front','callback','Addshapes');

function l=get_objects
l=get(gcbo,'userdata');
if strcmp(lower(get(l,'facecolor')),'none')
    x=get(l,'userdata');
    l=[l;x(:)];
end
    