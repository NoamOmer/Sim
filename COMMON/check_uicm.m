function check_uicm
a=get(gcbo,'parent');
set(get(a,'child'),'checked','off');
set(gcbo,'checked','on');
