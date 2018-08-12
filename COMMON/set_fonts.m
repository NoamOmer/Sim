function set_fonts(fig,axes_font_size,text_size)
if nargin==1
    axes_font_size=fig;
    text_size=fig;
    fig=gcf;
end
if nargin==2
    if ishandle(fig)
        text_size=axes_font_size;
    else
        text_size=axes_font_size; 
        axes_font_size=fig;
        fig=gcf;
    end
end
ax=findall(fig,'type','axes');
set(ax,'fontname','Times New Roman','fontsize',axes_font_size);
text_obj=findall(fig,'type','text');
set(text_obj,'fontname','Times New Roman','fontsize',text_size);
ax=findall(fig,'type','axes','tag','legend');
text_obj=findall(ax,'type','text');
set(text_obj,'fontname','Times New Roman','fontsize',8);
line_obj=findall(ax,'type','line');
set(line_obj,'markersize',5);
% set(ax,'box','off','xcolor','w','ycolor','w','color','none')
