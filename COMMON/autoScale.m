function autoScale(opt,fig)
a=findobj(fig,'type','axes');
ca=get(fig,'currentaxes');
if isempty(a),return;end
switch opt
    case 'GcaXY'
        compAS(ca,'XY')
    case 'GcaX'
        compAS(ca,'X')
    case 'GcaY'
        compAS(ca,'Y')
    case 'XY'
        compAS(a,'XY')
    case 'XYX'
        compAS(a,'X')
        compAS(ca,'Y')
    case 'XX'
        compAS(a,'X')
    case 'YY'
        compAS(a,'Y')
end
refresh(fig);

function compAS(AX,opt)
for o=1:length(AX)
    ax=AX(o);
    li=get(ax,'children');    
    if ~isempty(findstr(opt,'X'))
        xl=xlim;
        for k=1:length(li)
            switch lower(get(li(k),'type'))
            case {'line','patch'}
                x=get(li(k),'xdata');
                [Mi(k),Ma(k)]=min_max(x(:));
            otherwise                
                x=get(li(k),'position');
                Mi(k)=x(1);Ma(k)=x(1)+x(3);
            end
        end           
        xl=[min(Mi),max(Ma)];   
        if diff(xl)==0,xl=sort([xl(1)*.99,xl(2)*1.01]);end        
        set(ax,'xlim',xl)
    end
    if ~isempty(findstr(opt,'Y'))
        yl=ylim;
        for k=1:length(li)
            switch lower(get(li(k),'type'))
            case {'line','patch'}
                y=get(li(k),'ydata');
                [Mi(k),Ma(k)]=min_max(y(:));
            otherwise                
                y=get(li(k),'position');
                Mi(k)=y(2);Ma(k)=y(2)+y(4);
            end
        end           
        yl=[min(Mi),max(Ma)];           
        if diff(yl)==0,yl=sort([yl(1)*.99,yl(2)*1.01]);end    
        set(ax,'ylim',yl)
    end    
end