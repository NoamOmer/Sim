function setPointerInAxes(fig,pointerIn,pointerOut)
U=get(fig,'units');
set(fig,'units','normalized');
a=findobj(fig,'type','axes');
p=get(fig,'currentpoint');
x=0;
for k=1:length(a)
    P=get(a(k),'position');
    if p(1)>P(1) &p(1)<P(1)+P(3)
        if p(2)>P(2) &p(2)<P(2)+P(4)
            x=1;
        end
    end
end
if x==1
    set(gcf,'pointer',pointerIn);
else
    set(gcf,'pointer',pointerOut);
end
set(fig,'units',U);