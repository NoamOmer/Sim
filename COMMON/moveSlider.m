function moveSlider(slide)
if nargin==0
    slide=gcbo;
end
fig2=get(slide,'parent');
Sli=findobj(fig2,'tag','slider');
set(Sli,'enable','on');
prop=get(Sli,'userdata');
band=prop.displayInterval;
axe=findobj(fig2,'type','axes');

B=floor(get(Sli,'value'));
set(Sli,'value',B);
B=get(Sli,'value');
for k=1:length(axe)
    U=get(axe(k),'userdata');
    if length(U{1})>1
        Y=U{1};
        a=max([B,1]);
        a=min([a,length(Y)-band*U{3}]);
        b=a:min([length(Y),a+band*U{3}]);
        YY=Y(b);
        m=length(YY);
        set(U{4},'xdata',(b-1)/U{3},'ydata',YY),
        set(axe,'xlim',([b(1),b(end)]-1)/U{3})
    else
        [mes,errnum]=ferror(U{1});
        a=fseek(U{1},B*U{5},-1);
        if a<0
            if errnum==0
                [mes,errnum]=ferror(U{1});
                warndlg(mes);
            end
        else
            [YY,m]=fread(U{1},band*U{3},U{2});
            set(U{4},'xdata',(B+(1:m))/U{3},'ydata',YY),
            set(axe,'xlim',([1,m]+B)/U{3})
            fseek(U{1},0,-1);   
            ferror(U{1},'clear');
        end    
    end
end