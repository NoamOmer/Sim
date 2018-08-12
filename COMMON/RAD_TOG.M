function toggle_radio(tag,name)
if nargin==0
    b=gcbo;
elseif nargin==1
    b=tag;
else
    b=findobj(gcf,'tag',tag,'string',name);
end
tag=get(b,'tag');
name=get(b,'string');
col=get(b,'backgroundcolor');
a=findobj(gcf,'tag',tag);
set(a,'enable','on','value',0,'backgroundcolor',col);
if findstr(lower(get(b,'type')),'toggle')
    col=col*1.2;
end   
set(b,'value',1,'enable','inactive','backgroundcolor',col)

