function ax=get_axes(opt,handle)
switch opt
    case 'axes'
        fig=handle;
        ax=[];
        x=1;
        while ~isempty(x)
            [x,y,a]=get_points(fig,1);
            if ~isempty(a)
                ax(end+1)=a;
            end
        end
        return
    case {'InAxes','Outside','Right'}
        ax=handle;
        te='abcdefghij';        
        for k=1:length(ax)
            axUnits=get(ax(k),'units');
            set(ax(k),'units','normalized')
            p=get(ax(k),'pos');
            dy(k)=p(4);
            set(ax(k),'units',axUnits)
        end
        minDy=min(dy);
        maxDy=max(dy);
        factor=2*ones(size(ax));
        k=find(abs(dy-minDy)<100*eps);
        factor(k)=20;
        k=find(abs(dy-maxDy)<100*eps);
        factor(k)=8;
        dxFactor=10;
        if strmatch(opt,'Outside')
            dxFactor=-4;            
        end
        if strmatch(opt,'Right')
            dxFactor=.98;            
        end
        
        for k=1:length(ax)
            xl=get(ax(k),'xlim');
            x=xl(1)+diff(xl)/dxFactor;
            yl=get(ax(k),'ylim');
            y=yl(2)-diff(yl)/factor(k);
            t(k)=text(x,y,te(k),'parent',ax(k),'tag','index');
        end
        set(t,'fontw','bold')
        ax=t;
    case 'CopyLocation'
        ax=handle(:,1);
        destAx=handle(:,2);
        fig(1)=get(ax(1),'parent');
        fig(2)=get(destAx(1),'parent');
        figUnits{1}=get(fig(1),'units');
        figUnits{2}=get(fig(2),'units');        
        set(fig,'units','normalized');
        p=get(fig(1),'pos');
        set(fig(2),'pos',p);
        set(fig(1),'units',figUnits{1});      
        set(fig(2),'units',figUnits{2});
        
        for k=1:length(ax)
            axUnits{1}=get(ax(k),'units');
            set(ax(k),'units','normalized')
            axUnits{2}=get(destAx(k),'units');
            set(destAx(k),'units','normalized')
            
            p=get(ax(k),'pos');
            set(destAx(k),'pos',p)
            
            set(ax(k),'units',axUnits{1})
            set(destAx(k),'units',axUnits{2})
        end
        
end
