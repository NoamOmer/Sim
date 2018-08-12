function slider_band(opt,pointer)
% This function adds a slider under an axis, which allows you to browse through the plot.
% By changing an editable box (below the slider) you may change the Xzoom property of the axes 
% (the range of the x limit). the initial range is 10% of the lines ranges. For example, if the lines 
% within the axes have x ranges of: [0,100] and [-100,50] , the total range is 200, thus the 
% initial range will be 20.
% Calling the function with no arguments will draw the slider below the current axis of the current figure.
% The first argument, opt, should be a string of 'create' or an empty string.
% The second argument, pointer, may be either an handle to a figure or an array of handles to axes.
% Note: if pointer is an array of axes the slider will act simultaneously on all of them.
%     slider_band;
%     slider_band('create')
%     slider_band('create',fig)
%     slider_band('create',[axes1,axes2,...])
    
if nargin==0 | isempty(opt)
    opt='create';
end
if nargin<2
    fig=gcbf;
else
    fig=pointer;
end
if isempty(fig)
    fig=gcf;
end
if strcmp(get(fig(1),'type'),'axes')
    ax=fig;
    fig=get(ax(1),'parent');
else
    ax=get(fig,'currentaxes');
end
for k=1:length(ax)
    p(k,:)=get(ax(k),'pos');% collecting the positions of all axes
end
[a,I]=sort(p(:,2));% sorting the positions from lowest to highest axes
ax=ax(I);
l=findobj(ax(1),'type','line');% all lines within the lowest axes
for k=1:length(l)
    x{k}=get(l(k),'xdata');
    y{k}=get(l(k),'ydata');
    [miX(k),maX(k)]=min_max(x{k});    
    [miY(k),maY(k)]=min_max(y{k});    
end
minX=min(miX);
maxX=max(maX);
switch opt
    case 'create'
        Child=get(fig,'children');
        CurUn=get(Child,'units');% storing the units of all objects
        set(Child,'units','cent');% changing the units to cent
        set(fig,'units','cent')
        p=get(fig,'pos');
        p(2)=p(2)-1;p(4)=p(4)+1;
        set(fig,'pos',p)% we increased the figure height by 1 cm
        set(fig,'units','norm')
        for k=1:length(Child)
            p=get(Child(k),'pos');
            p(2)=p(2)+1;
            set(Child(k),'pos',p);% we increased the bottom posision of all objects by 1 cm
            % restoring the units
            if iscell(CurUn)
                set(Child(k),'units',CurUn{k});
            else
                set(Child(k),'units',CurUn);
            end
        end
        axUnits=get(ax(1),'units');
        set(ax(1),'units','cent');
        pCent=get(ax(1),'pos');
        set(ax(1),'units','norm');
        pNorm=get(ax(1),'pos');
        set(ax(1),'units',axUnits)
        factor=pNorm(2)/pCent(2);
        p=pNorm;
        y=p(2)-factor*[2,1.5];
        sliderPos=[p(1),y(2),p(3),.03];
        % drawing the slider
        slider=uicontrol('Style','slider','units','norm','min',minX,'max',minX+(maxX-minX)*.95,...
            'position',sliderPos,'sliderstep',[.1 .9]/5,'tag','slider','callback','slider_band(''slider'');',...
            'value',minX,'userdata',ax);        
        % drawing the range edit box
        h(1)=add_button('text',[p(1)+p(3)/3,y(1),.05,.04],'<---',[],'slider_text');
        h(2)=add_button('edit',[p(1)+p(3)/3+.05,y(1),p(3)/3-.1,.04],num2str((maxX-minX)/10,2),...
            [1 1 1],'band','slider_band(''band'');');
        h(3)=add_button('text',[p(1)+2*p(3)/3-.05,y(1),.05,.04],'--->',[],'slider_text');
        h(4)=add_button('check',[p(1),y(1),p(3)/3-.02,.04],'Autoscale Y',[],'autoscale',...
            'slider_band(''autoscale'');');
        set(h,'fontsize',8);
        opt='band';
end
BAND=findobj(fig,'tag','band');
band=str2num(get(BAND,'string'));
Slider=findobj(fig,'tag','slider');
set(Slider,'enable','on');
ax=get(Slider,'userdata');
switch opt
    case 'band'
        if band<0
            band=maxX-minX;
        end
        if maxX-minX<=band
            set(Slider,'enable','inactive');
            set(ax,'xlim',[minX,maxX])
            set(BAND,'string',num2str(maxX-minX,2))
        end
        if strcmp(lower(get(Slider,'enable')),'on')
            xl=get(ax(1),'xlim');
            set(ax,'xlim',xl(1)+[0,band])
            set(Slider,'min',minX,'max',maxX-band,'sliderstep',[.1,.9]*band/(maxX-minX),...
                'value',max([minX,xl(1)]))
        end
    case 'slider'
        B=get(Slider,'value');
        set(ax,'xlim',B+[0,band])
    case 'reset_slider'
        set(ax,'xlim',minX+[0,band])
        slider_band('band',fig);       
    case 'update'
        set(ax,'xlim',minX+[0,band])    
        slider_band('band',fig); 
    case 'delete'
        delete(Slider);
        delete(BAND);
        delete(findobj(fig,'tag','slider_text'));
        delete(findobj(fig,'tag','autoscale'));
        return
end
% checking the autoscale flag
AS=get(findobj(fig,'tag','autoscale'),'value');
for k=1:length(ax)
    set_scale(ax(k),AS,miY,maY);
end

function  set_scale(ax,AS,miY,maY)
l=findobj(ax,'type','line');

xl=get(ax,'xlim');
miY=[];
maY=[];
for k=1:length(l)
    x=get(l(k),'xdata');
    y=get(l(k),'ydata');
    if AS==1
        L=find(x>=xl(1) & x<=xl(2));
        if ~isempty(L)
            [miY(end+1),maY(end+1)]=min_max(y(L));            
        end
    else
        [miY(end+1),maY(end+1)]=min_max(y);
    end
end
minY=min(miY);
maxY=max(maY);
set(ax,'ylim',[minY,maxY])