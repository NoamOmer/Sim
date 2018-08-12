function x=createWaitBar(opt,val,WB);
if nargin==3
    fig=WB.figure;
    figure(fig);
    x=get(fig,'userdata');
end
switch opt
    case 'init'
        a=findall(0,'tag','WaitBar');
        if ~isempty(a)
            delete(a)
        end
        x=[];
        if nargin==1
            val=0;
        end      
        x.figure=figure('visible','off','units','normalized','pos',[.3,0.45,0.4,0.1],'number','off','menu','none');
        set(x.figure,'units','cent');
        pos=get(x.figure,'pos');
        pos(4)=1+1*length(val);
        set(x.figure,'pos',pos);        
        wb=x.figure;
        set(wb,'handlevisibility','on','closerequestfcn','closereq','tag','WaitBar')
        set(0,'currentfigure',wb)
        L=length(val)+1;
        for k=1:length(val)
            x.axis(k)=axes('units','normalized','pos',[0.05,(L-k)/pos(4),0.9,.5/pos(4)]);
            x.title(k)=get(x.axis(k),'title');
            x.patch(k)=patch([0 0 val(k) val(k) 0]/100,[0 1 1 0 0],'r');
            x.text(k)=text(0.45,0.5,'0 %','parent',x.axis(k));
        end
        set([x.text,x.title],'fontweight','bold','color',[0 0 1],'Interpreter','none')
        set(x.axis,'drawmode','fast','box','on','xlim',[0 1],'ylim',[0 1],'xtick',[],'ytick',[])
        set(x.figure,'userdata',x,'visible','on')        
        add_button('push',[0.3,0.01,0.4,.6/pos(4)],'Press to terminate',[],'','createWaitBar(''terminate'')');
        alignTitles(x.figure,[0.5,.85,0])
        x.error = 0;
        drawnow
    case 'terminate'
        x=get(gcbf,'userdata');
        delete(x.figure)
		x.error = 1;
        error('Terminated by user');
        return
    case 'percent'
        for k=1:length(val)        
            set(x.text(k),'string',[num2str(round(val(k))),'%'])
            set(x.patch(k),'xdata',[0 0 val(k)/100 val(k)/100 0])        
        end
    case 'Addpercent'             
        for k=1:length(val)
            X=max(get(x.patch(k),'xdata')); 
            set(x.text(k),'string',[num2str(round((X+val(k)/100)*100)),'%'])
            set(x.patch(k),'xdata',[0 0 X+val(k)/100 X+val(k)/100 0])        
        end
    case 'title'
        for k=1:length(x.axis)
            if iscell(val)
                set(x.title(k),'string',val{k})
            else
                set(x.title(k),'string',val)
            end
        end
end
drawnow;
% figure(x.figure);
% set(0,'currentfigure',x.figure)