function toggle_children(obj);
if nargin<1
    obj=gcbo;
end
child=get(obj,'userdata');
if get(obj,'value')==1
    set(child,'enable','on');
else
    set(child,'enable','inactive')
end