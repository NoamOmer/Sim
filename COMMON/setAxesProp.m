function setAxesProp(fig,LabelFontName,LabelFontSize,AxesFontsize)

a=findobj(fig,'type','axes');
for l=1:length(a)
    s=get(a(l),'xlabel');
    set(s,'fontname',LabelFontName,'fontsize',LabelFontSize)
    s=get(a(l),'ylabel');
    set(s,'fontname',LabelFontName,'fontsize',LabelFontSize)
    s=get(a(l),'title');
    set(s,'fontname',LabelFontName,'fontsize',LabelFontSize)
end
