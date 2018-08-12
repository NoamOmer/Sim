function new_fig=export_fig(fig,settings_file)
if nargin<1
    fig=gcf;
end

if nargin<2
    settings_file=get_settings_file;
end
if isempty(settings_file)
    return
elseif isstruct(settings_file)
    S=settings_file;
elseif exist(settings_file)
    load (settings_file,'-mat');
else 
    return
end

fig_type=get_fig_type(fig);
switch fig_type
    case '3D plot'
        newfig=figure;
        set(newfig,'units','cent','position',S.paperpos,'paperuni','cent','paperpos',S.paperpos);
        colormap(S.map);
        AxesTag={'3dplot';'2dplot';'hbar'};
        %new position for 3 axes
        s1=.7;
        x0=.11+[0,0,s1+.02];
        y0=.13+[0,s1+.02,0];
        dx=[s1,s1,.03];
        dy=[s1,.98-s1-.15,s1];
        for k=1:3
            ax(k)=findobj(fig,'type','axes','tag',AxesTag{k});
            NewAxes(k)=copyobj(ax(k),newfig);
            set(NewAxes(k),'units','normalized','position',[x0(k),y0(k),dx(k),dy(k)],...
                'xcolor',S.axeslinecolor,'ycolor',S.axeslinecolor,'box','on','layer','top')
            xl(k,1:2)=get(NewAxes(k),'xlim');
            yl(k,1:2)=get(NewAxes(k),'ylim');
            x_label(k)=get(NewAxes(k),'xlabel');
            y_label(k)=get(NewAxes(k),'ylabel');
        end
        set([NewAxes,x_label,y_label],...
            'color',S.textcolor,'fontsize',S.axesfontsize,'fontname','Times New Roman')

        set(NewAxes(2),'xlim',xl(1,:),'xtick',[],'color',S.plot2dbackcolor);
        vc=findobj(NewAxes(2),'type','line');
        if ~isempty(vc)
            set(vc,'color',S.plot2dlinecolor,'ButtonDownFcn','');
        end
        if strcmp(lower(get(fig,'userdata')),'hosa')
            xlabelStr='Frequency [Hz]';
        else
            xlabelStr='Time [sec]';
        end
        set(y_label(1),'string','Frequency [Hz]','color',S.textcolor)
        set(x_label(1),'string',xlabelStr,'color',S.textcolor)
        set(y_label(2),'string','BPM^2/Hz')

        pa=findobj(newfig,'type','patch');
        set([NewAxes pa(:)'],'ButtonDownFcn','');
        vc=get(pa,'cdata');
        v=cat(1,vc{:});
        map=get(newfig,'colormap');
        lm=size(map,1);
        co=round((v-min(v))/(max(v)-min(v))*lm);
        newco=min(co+3,lm);
        for k=1:length(pa)
            set(pa(k),'edgecolor',map(newco(k),:));
        end
        set(newfig,'menu','figure')
        drawnow
        change_labels(NewAxes);        
    otherwise
        newfig=figure;
        set(newfig,'units','cent','position',S.paperpos,'paperunits','cent','paperpos',S.paperpos)
        colormap(S.map);
        ax=findobj(fig,'type','axes');   
        N_axes=length(ax);
        for k=1:N_axes
            tt=get(ax(k),'position');
            t(k)=tt(2);
            naxe(k)=copyobj(ax(k),newfig);
        end
        [N,I]=sort(t);I=I(end:-1:1);
        ax=ax(I);naxe=naxe(I);
        s1=0.75;
        step=.75/N_axes;
        ypos=.1+step*((1:N_axes)-1)*1.1;
        ypos=ypos(end:-1:1);
        for k=1:length(ax)        
            set(naxe(k),'pos',[.11,ypos(k),s1,step],'xcolor',S.axeslinecolor,'ycolor',...
                S.axeslinecolor,'box','on','layer','top')
        end
        pa=findobj(naxe,'type','patch');
        set([naxe pa(:)'],'ButtonDownFcn','');
        map=get(newfig,'colormap');
        set(findobj(naxe,'type','line'),'color',S.plot2dlinecolor)
        for k=1:length(ax)
            a(1)=get(naxe(k),'xlabel');
            a(2)=get(naxe(k),'ylabel');
            a(3)=get(naxe(k),'title');
            set(a,'fontsize',S.textfontsize,'fontname','timesnewroman','color',S.textcolor)
        end
        set(naxe,'fontsize',S.axesfontsize,'color',S.plot2dbackcolor)
        set(fig,'menu','figure')
        drawnow
        change_labels(naxe);
        a=findall(newfig,'type','line');
        set(a,'buttondownfcn','');
        a=findall(newfig,'type','axes');
        set(a,'buttondownfcn','');        
end



function fig_type=get_fig_type(fig)
fig_type=get(fig,'tag');

function change_labels(ax)
a=questdlg('Would you like to change labels and titles ?','','Yes','No','Yes');
if strmatch(a,'No')
    return
end
fig=figure('units','normalized','position',[.2,.2,.4,.4],'windowstyle','modal');
tex=add_button('text',[.01,.8,.98,.1]);
a(1)=add_button('text',[.05,.6,.4,.1],'Title :');
a(2)=add_button('text',[.05,.48,.4,.1],'Xlabel :');
a(3)=add_button('text',[.05,.36,.4,.1],'Ylabel :');

b(1)=add_button('edit',[.5,.6,.4,.1]);
b(2)=add_button('edit',[.5,.48,.4,.1]);
b(3)=add_button('edit',[.5,.36,.4,.1]);

c(1)=add_button('push',[.2,.1,.25,.1],'OK',[],'','set(gcf,''userdata'',1);');
c(2)=add_button('push',[.55,.1,.25,.1],'Cancel',[],'','set(gcf,''userdata'',0);');
%ax=ax(end:-1:1);
for k=1:length(ax)
    t=get(ax(k),'title');
    x=get(ax(k),'xlabel');
    y=get(ax(k),'ylabel');
    set(tex,'string',['Enter labels for axis number ',num2str(k)]);
    set(b(1),'string',get(t,'string'));
    set(b(2),'string',get(x,'string'));
    set(b(3),'string',get(y,'string'));
    
    waitfor(fig,'userdata')
    if get(fig,'userdata')==0
        delete(fig)
        return
    end
    set(t,'string',get(b(1),'string'));
    set(x,'string',get(b(2),'string'));
    set(y,'string',get(b(3),'string'));
    set(fig,'userdata',[])
end    
delete(fig)

function settings_file=get_settings_file
a1=which('spat');
d=[get_file_prop(a1,'path'),'export'];
A=dir([d,'\*.set']);
a=[];
if ~isempty(A)
    for k=1:length(A)
        S{k}=A(k).name;
    end
    S{k+1}='Change parameters';
    S{k+2}='Create a new set';
    FigureUnits=get(0,'defaultfigureunits');
    set(0,'defaultfigureunits','pixels');
    a=listdialog('liststring',S,'name','Export set','selectionmode','single');
    set(0,'defaultfigureunits',FigureUnits);
end
if isempty(a)
    return
elseif a==k+1;
    [f,p]=uigetfile('*.set');
    if p==0,return,end
    settings_file=createNewSet([p,f]);        
elseif a==k+2 
    settings_file=createNewSet;        
else
    settings_file=[d,'\',S{a}];
end
