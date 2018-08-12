function fig2=plot_3d(varargin)
%
% This function plots a 2d matrix (used for time-freq analysis and HOSA)
% If P is complex, it plots the matrix B=abs(P) (the squared absolute value)!!!!
%
% If time series is given a third axes is open for it.
%
% Pressing on the main axes gives the t f A values.
% A line is drawn in the time series axes at the same t.
% For a squared matrix (HOSA) another 2 line are drawn at f and t+f.    
%
% Rectangles/Polygons: 
%                 (x,y) refares to the mean x and y.
%                 Max point is the point with maximum value inside the curve and its phase.
%                 Mean value is the integral inside the curve P(x,y)*dx*dy normalized by the
%                 number of points inside the polygon.
%
% Pressing on a curve will highlight the corresponding time interval 
% of the time series. 
%
% The Shading checkbox is effective only for Surf, Mesh and Waterfall.
% The Levels button is effective only for Contour and Contourf.
% The Threshold level is used only for display.

K=varargin;
if length(K)>0
    if length(K{1}(:))==1
        fig2=K{1};
        if ~ishandle(fig2) | isempty(strcmp(get(fig2,'tag'),'3D plot'))
            fig2=gcf;
        end
        choice=K{2};
    else
        fig2=[];
    end
end

series=[];lx=[];ly=[];INFOTXT='';figureName='3D plot';ff=[];tt=[];FileN=[];
if isempty(fig2)
    ff=[];tt=[];l=[];
    for k=1:length(K)
        if isstr(K{k})
            if findstr(lower(K{k}),'info')
                INFOTXT=K{k+1};
                l=k;
            end
        end
    end
    if ~isempty(l),    K(l:l+1)=[];end
    l=[];
    for k=1:length(K)
        if isstr(K{k})
            if findstr(lower(K{k}),'filename')
                FileN=K{k+1};
                l=k;
            end
        end
    end
    if ~isempty(l),    K(l:l+1)=[];end
    
    P=K{1};choice='Refresh';
    K(1)=[];
    if ~isempty(K)
        if ~isstr(K{1}),ff=K{1};K(1)=[];end
    end
    if ~isempty(K)
        if ~isstr(K{1}),tt=K{1};K(1)=[];end
    end
    if ~isempty(K)
        if ~isstr(K{1}),
            series=K{1};K(1)=[];
            [m,n]=size(series);
            time_series=1;
            if n==1 |n==3,
                series(:,end+1)=(((1:length(series(:,1)))-1)*(tt(2)-tt(1)))';
            end
        end
    end
    if ~isempty(K)
        figureName=K{1};
    end
    
    [m,n]=size(P);
    if isempty(ff),ff=1:m;end
    if isempty(tt),tt=1:n;end   
    if isempty(series),time_series=0;end
    fig2=add_figure(time_series,figureName,INFOTXT,FileN);
    set(fig2,'userdata','TFD')
    set(findobj(fig2,'tag','add_point'),'enable','off')
    if size(ff)==size(tt)
        if ff==tt
            set(fig2,'userdata','HOSA')
            set(findobj(fig2,'tag','add_point'),'enable','on')
        end
    end

    CLOSE=findobj(fig2,'tag','close');
    [m,n]=size(P);
    if length(ff)==n &length(tt)==m
        P=P';
    end
    ly=1:length(ff);
    lx=1:length(tt);
    lx=lx(:);ly=ly(:);ff=ff(:);tt=tt(:);
    set(CLOSE,'userdata',[{P},{ff},{tt},{series}]);
    set(findobj(fig2,'tag','3dplot'),'xlim',[tt(1),tt(end)],'ylim',[ff(1),ff(end)])
    if ~isempty(findstr(lower(figureName),'wavelet')), 
        set(findobj(fig2,'tag','drawMode'),'value',1);
        set(findobj(fig2,'tag','fFactor'),'value',1);
    end
    if ~isempty(findstr(lower(figureName),'bicoherence')),
        set(findobj(fig2,'tag','drawMode'),'value',3);
    end
    if ~isempty(findstr(lower(figureName),'skewness')), 
        set(findobj(fig2,'tag','drawMode'),'value',3);
    end        
end
CLOSE=findobj(fig2,'tag','close');
add_rect=findobj(fig2,'tag','add_rect');
add_poly=findobj(fig2,'tag','add_poly');
add_text=findobj(fig2,'tag','add_text');
info=findobj(fig2,'tag','info');
p3d=findobj(fig2,'tag','3dplot');
p2d=findobj(fig2,'tag','2dplot');


set(fig2,'pointer','watch');
switch choice
case 'Refresh'
    plot_fcn(fig2,add_rect,add_poly,add_text);
case 'chaxes'% Change Frequency limits
    a=gcbo;
    [Fx,Fy]=SelectAxesLimits(p3d);
    set(a,'userdata',[Fx(:)';Fy(:)']);
case 'hbarlevels' % Change zlim of the main axes
    hbar=findobj(fig2,'tag','hbar');
    pp=get(p3d,'userdata');    
    minP=pp(1);
    maxP=pp(2);
    levels=str2num(get(findobj(fig2,'tag','levels'),'string'));
    hlevels=selectZlim(minP,maxP);
    set(hbar,'userdata',hlevels);
case {'arect','srect'} % add rectangle
    if strcmp(choice,'arect')
        r=getrect(fig2);
        pos.x=r(1)+[0,r(3),r(3),0,0];
        pos.y=r(2)+[0,0,r(4),r(4),0];
    elseif strcmp(choice,'srect')
        [x,y]=SelectRectLimits;
        if isempty(x),
            return;
        end
        pos.x=[x(1),x(2),x(2),x(1),x(1)];
        pos.y=[y(1),y(1),y(2),y(2),y(1)];
    end   
    [MAX,Phase,VOL,peak]=comp_rectangle(fig2,pos);
    RECT=get(add_rect,'userdata');
    RECT.POS{end+1}=pos;
    RECT.PEAK{end+1}=peak;
    RECT.MAX{end+1}=MAX;
    RECT.Phase{end+1}=Phase;
    RECT.VOL{end+1}=VOL;
    set(add_rect,'userdata',RECT);
    draw_line(fig2,add_rect,'temp_rect'); 
case 'poly' %add polygon
    [pos.x,pos.y]=getline(fig2);
    pos.x(end)=pos.x(1);pos.y(end)=pos.y(1);     
    [MAX,Phase,VOL,peak]=comp_polygon(fig2,pos);
    RECT=get(add_poly,'userdata');
    RECT.POS{end+1}=pos;
    RECT.PEAK{end+1}=peak;
    RECT.MAX{end+1}=MAX;
    RECT.Phase{end+1}=Phase;
    RECT.VOL{end+1}=VOL;
    set(add_poly,'userdata',RECT);
    draw_line(fig2,add_poly,'temp_poly');       
case {'wavInt','awavInt'} % Select regions for wavelet integrals
    if strmatch(choice,'awavInt')
        [LF,HF]=addWavIntFig;
    else
        [x,y]=get_points(p3d,4);
        if length(y)==4
            y=sort(y);
            LF=y(1:2);HF=y(3:4);
        end
    end
    if length(LF)==2 & length(HF)==2
        U=get(CLOSE,'userdata');
        P=U{1};ff=U{2};tt=U{3};series=U{4};
        X=get(p3d,'xlim');
        K=find(tt>X(1) &tt<X(2));
        K1=find(series(:,2)>X(1) &series(:,2)<X(2));
        show_integrals(series(K1,2),series(K1,1),tt(K),ff,P(:,K),LF,HF,get(findobj(gcbf,'tag','INF'),'userdata'),...
            get(findobj(gcbf,'tag','FileName'),'string'));
    end
case 'clear' % clear markers 
    clearBut=findobj(fig2,'tag','clear');
    opt=get(clearBut,'userdata');
    POLY=struct('POS',[],'PEAK',[],'MAX',[],'Phase',[],'VOL',[]);
    P=struct('POS',[],'TEXT',[]);
    if opt==1
        delete(findobj(fig2,'tag','temp_line'))   
        delete(findobj(fig2,'tag','temp_rect'))   
        delete(findobj(fig2,'tag','temp_poly'))   
        delete(findobj(fig2,'tag','temp_point'))   
        delete(findobj(fig2,'tag','temp_text'))   
        set(add_poly,'userdata',POLY);
        set(add_rect,'userdata',POLY);        
        set(add_text,'userdata',P);        
    elseif opt==2
        set(add_rect,'userdata',POLY);
        delete(findobj(p3d,'type','line','tag','temp_rect'));   
        delete(findobj(fig2,'tag','temp_line'))   
    elseif opt==3
        set(add_poly,'userdata',POLY);
        delete(findobj(p3d,'type','line','tag','temp_poly'));   
        delete(findobj(fig2,'tag','temp_line'))   
    elseif opt==4
        delete(findobj(fig2,'tag','temp_point'))
    elseif opt==5
        delete(findobj(fig2,'tag','temp_text'));
        set(add_text,'userdata',P);        
    end
case 'inf' %display object information
    line_inf=get(gco,'userdata');
    x=get(gco,'xdata');
    y=get(gco,'ydata');
    MAX=line_inf(1);
    Phase=line_inf(2);
    VOL=line_inf(3);
    px=line_inf(4);   py=line_inf(5);
    string=['(',num2str(px,3),',',num2str(py,3),'):   MAX point: value  ',num2str(MAX,3),...
            ',   Phase  ',num2str(Phase,3),',   Volume: ',num2str(VOL,3)];
    set(info,'string',string)   
    if ~isempty(p2d)
        p2d_line=findobj(p2d,'tag','xy_2dplot');
        delete(findobj(fig2,'tag','temp'));
        X=get(p2d_line,'xdata');
        Y=get(p2d_line,'ydata');
        if iscell(X),
            X1=X{1};Y1=Y{1};X2=X{2};Y2=Y{2};p2d_line=p2d_line(1);
        else
            X1=X;X2=X;Y1=Y;Y2=Y;
        end
        lx=find(X1>min(x) &X1<max(x));
        ly=find(X2>min(y) &X2<max(y));
        set(fig2,'currentaxes',p2d);
        line(X1(lx),Y1(lx),'color',[1 1 1],'tag','temp','linewidth',2);
        set(fig2,'currentaxes',p3d);
    end
case 'location' %add arrow on the current point
    p=get(gca,'currentpoint');
    add_current_location(fig2,p(1,1),p(1,2),info);
case 'text' %add text
    [a,b]=get_points(fig2,1);
    te=inputdlg('Enter text in the edit box');
    c=text(a,b,te);
    P=get(add_text,'userdata');
    P.POS{end+1}=[a,b];
    P.TEXT{end+1}=te;
    set(add_text,'userdata',P)
    if get(fig2,'userdata')==1,   col=[0 0 0];else,col=[1 0 0];end
    set(c,'color',col,'tag','temp_text');
case 'qpc'
    if strcmp(get(gcbo,'type'),'axes')
        a=get(p3d,'currentpoint');
        a=[a(1,1),a(1,2)];
    else
        a=get_points(p2d,2);
        if max(size(a))<2,        return;end
        a=a(:,1);
    end
    add_QPC(fig2,a(1),a(2),info);
case 'information'
    h=msgbox(get(gcbo,'userdata'),'createMode','modal'); 
case 'zoomIn'
    r=getrect(fig2);
    if r(3)>0 & r(4)>0
        x=get(p3d,'xlim');
        y=get(p3d,'ylim');
        zi=findobj(fig2,'tag','zoomIn');
        zo=findobj(fig2,'tag','zoomOut');
        U=get(zi,'userdata');
        U{end+1}=[x;y];
        set(zi,'userdata',U)
        set(zo,'enable','on')
        x=r(1)+[0,r(3)];
        y=r(2)+[0,r(4)];
        set(p3d,'xlim',x,'ylim',y);
        set(p2d,'xlim',x);
    end
case 'zoomOut'
    zi=findobj(fig2,'tag','zoomIn');
    zo=findobj(fig2,'tag','zoomOut');
    U=get(zi,'userdata');
    X=U{end};
    set(p3d,'xlim',X(1,:),'ylim',X(2,:));
    set(p2d,'xlim',X(1,:));
    U(end)=[];
    set(zi,'userdata',U)
    if isempty(U)
        set(zo,'enable','off')
    end
case 'zoom restore'
    aLim=findobj(fig2,'tag','axesLimit');
    X=get(aLim,'userdata');
    xl=X(1,:);
    yl=X(2,:);   
    set(p3d,'xlim',xl,'ylim',yl);
    set(p2d,'xlim',xl);    
case 'Filter'
  a1=inputdlg({'Enter Filter size. Use m,n format.'},'a',1,{'25,1'});
  a=num2str(['[',a1{1},']']);
   k=get(gco,'uicontextmenu');
   k1=get(k,'children');
   set(k1,'checked','off')
   l=uimenu(k,'label',a1{1},'callback',['check_uicm;set(gco,''userdata'',',a,');']...
       ,'checked','on');
   k1=[k1(1);l;k1(2:end)];
   set(k,'children',k1);
   set(gco,'userdata',str2num(a));
case 'save'
    save_data(fig2);
case 'exit'
    exit_quest(fig2);
end
%if ishandle(fig2),refresh(fig2);end
set(findobj(0,'type','figure'),'pointer','arrow');
return
%************************************************************************************   
function plot_fcn(fig2,add_rect,add_poly,add_text);

CLOSE=findobj(fig2,'tag','close');
U=get(CLOSE,'userdata');
P=U{1};ff=U{2};tt=U{3};series=U{4};
p3d=findobj(fig2,'type','axes','tag','3dplot');
aLim=findobj(fig2,'tag','axesLimit');
X=get(aLim,'userdata');
if isempty(X)
    xl=min_max(tt);
    yl=min_max(ff);
    set(aLim,'userdata',[xl(:)';yl(:)'])
else
    xl=X(1,:);
    yl=X(2,:);   
end
lx=find(tt>=xl(1) & tt<=xl(2));
ly=find(ff>=yl(1) & ff<=yl(2));
P=P(ly,lx);tt=tt(lx);ff=ff(ly);
tt=tt(:);ff=ff(:);
BC=abs(P);

S=['plot_3d(',num2str(fig2),','];
t_level=str2num(get(findobj(fig2,'tag','t_lev'),'string'));
l=find(BC<t_level);BC(l)=0;
if get(findobj(fig2,'tag','log'),'value')==1
    kk=find(BC<=0);
    k=find(BC>0);
    TT=get(findobj(fig2,'tag','log'),'userdata');
    if findstr(TT,'log2'),   BC(k)=log2(BC(k));
    elseif findstr(TT,'log10'), BC(k)=log10(BC(k));
    else BC(k)=log(BC(k));end           
    mm=min(BC(:));
    BC(kk)=mm;
end
if get(findobj(fig2,'tag','fFactor'),'value')==1
    TT=get(findobj(fig2,'tag','fFactor'),'userdata');
    switch TT
    case 'f',  BC=BC.*(ff*ones(size(tt')));
    case 'f2',  BC=BC.*((ff.*ff)*ones(size(tt')));
    case 'f3',  BC=BC.*(sqrt(ff)*ones(size(tt')));
    end
end

levels=str2num(get(findobj(fig2,'tag','levels'),'string'));
p3d=findobj(fig2,'tag','3dplot');
p2d=findobj(fig2,'tag','2dplot');
hbar=findobj(fig2,'tag','hbar');
add_point=findobj(fig2,'tag','add_point');
clear_markers=findobj(fig2,'tag','clear');

plot_method=get(findobj(fig2,'tag','drawMode'),'value');
if ~isempty(p2d)
    set(fig2,'currentaxes',p2d);
    delete(findobj(p2d,'type','line'));
    col=[.9 .9 0;.9 .3 0;0 .9 0];
    [m,n]=size(series);
    serTime=series(:,n);
    l=find(serTime>tt(1) &serTime<tt(end));
    for i=1:n-1
        line(serTime(l),series(l,i),'tag','xy_2dplot','color',col(i,:),'linewidth',1.5);
    end
    axis tight
end   
set(fig2,'currentaxes',p3d)
sh=get(findobj(fig2,'string','Shading'),'value');
hbarLevels=get(hbar,'userdata');
if ~isempty(hbarLevels)
    mm=min_max(BC(:));
    if hbarLevels(1)<mm(2) & hbarLevels(2)>mm(1)
        levels=linspace(hbarLevels(1),hbarLevels(2),levels);
    end
end
fil=findobj(fig2,'tag','filter');
if get(fil,'value')==1
    len=get(fil,'userdata');
    l=len(1);
    Bfil=2.^(0:(l-1)/2);
    Bfil=[Bfil,2^((l+1)/2),Bfil(end:-1:1)];
    if len(2)==1
        BC=cmed(BC',l)';
    else
        Bfil=Bfil'*Bfil;
        Bfil=Bfil./sum(Bfil(:));
        BC=filter2(Bfil,BC);
    end
end
        
switch plot_method
case 1
    if length(levels)==1
        levels=min_max(BC(:));
    end    
    imagesc(tt,ff,BC,min_max(levels)'),axis xy
case 2
    contour(tt,ff,BC,levels);
case 3
    contourf(tt,ff,BC,levels);   
%    conf(tt,ff,BC,levels);   
case 4
    if length(levels)>1
        k=find(BC(:)>levels(end));
        BC(k)=levels(end);
        k=find(BC(:)<levels(1));
        BC(k)=levels(1);
    end
    surf(tt,ff,BC);
case 5
    mesh(tt,ff,BC);
case 6
    waterfall(tt,ff,BC);
end

minP=min(min(BC));
maxP=max(max(BC));
K=findobj(fig2,'tag','method');
set(K,'backgroundcolor',[.1529,.0824,.5412])
for k=1:length(K)
    if ~isempty(findstr(get(K(k),'callback'),plot_method))
        set(K(k),'backgroundcolor',[.4,.2,1])
    end
end    
if ~isempty(find([4,5,6]==plot_method))
    if sh==1,shading interp;end
    set([add_rect(:);add_poly(:)],'enable','off')
    set([add_point(:);clear_markers(:)],'enable','off')
    set(p3d,'tag','3dplot','drawmode','fast','userdata',[minP,maxP]);
    colorbar(hbar);
    set(hbar,'tag','hbar','userdata',hbarLevels,'buttondownfcn',[S,'''hbarlevel'');']);
else
    set([add_rect(:);add_poly(:)],'enable','on')
    if strcmp(get(fig2,'userdata'),'HOSA')
        set([add_point(:);clear_markers(:)],'enable','on')
    else
        set(clear_markers(:),'enable','on')
    end
    colorbar(hbar);
    set(hbar,'tag','hbar','userdata',hbarLevels,'buttondownfcn',[S,'''hbarlevels'');']);
    set(get(hbar,'children'),'buttondownfcn',[S,'''hbarlevels'');']);
    set(p3d,'tag','3dplot','drawmode','fast','buttondownfcn',[S,'''qpc'');'],...
        'userdata',[minP,maxP]);    
    s=findobj(p3d,'type','patch');
      set(s,'hittest','off');
%     set(s,'buttondownfcn',[S,'''location'');']);
    s=findobj(p3d,'type','image');
%     set(s,'buttondownfcn',[S,'''location'');']);
    set(s,'hittest','off');
    set(p2d,'buttondownfcn',[S,'''location'');']);
    s=findobj(p2d,'type','line');
    set(s,'buttondownfcn',[S,'''location'');']);
    draw_line(fig2,add_rect,'temp_rect'); 
    draw_line(fig2,add_poly,'temp_poly');   
    add_prev_text(fig2,add_text);
end

%************************************************************************************   
%************************************ The main figure
%************************************************************************************   
function fig2=add_figure(time_series,figureName,INFOTXT,FileN)
fig2=figure('name',figureName,'tag','3D plot','menu','figure','number','off','units','normalized',...
    'position',[.1,.06,.8,.9],'color',.75*[1,1,1],'renderer','zbuffer');
POLY=struct('POS',[],'PEAK',[],'MAX',[],'Phase',[],'VOL',[]);
P=struct('POS',[],'TEXT',[]);
l=.83+[0,.01,.05];
a(1)=add_button('frame',[.01,l(1),.17,.04*2+.015]);
S=['plot_3d(',num2str(fig2),','];
add_button('text',[.02,l(3),.15,.035],'Draw method');
add_button('popup',[.02,l(2),.15,.035],{'Image','Contour','Contourf','Surf','Mesh','Waterfall'},...
    [1 1 1],'drawMode','',[],'value',3);

umenu=uicontextmenu;
frac=uicontextmenu;
addrect=uicontextmenu;
clearMark=uicontextmenu;
wavIntegral=uicontextmenu;
zoom_restore=uicontextmenu;
filMenu=uicontextmenu;
l=.55+[0,.01,.05,.09,.13,.17,.21];
a(2)=add_button('frame',[.01,l(1),.17,.04*6+.015]);
add_button('checkbox',[.02,l(7),.15,.035],'Shading');
add_button('check',[.02,l(6),.15,.035],'log values',[],'log','','log',...
    'uicontextmenu',umenu);
add_button('check',[.02,l(5),.15,.035],'Freq. factor',[],'fFactor','','f',...
    'uicontextmenu',frac);
add_button('check',[.02,l(4),.15,.035],'Use filter',[0 1 1],'filter','',[7,7],...
    'foregroundcolor',[1 1 1],'uicontextmenu',filMenu);

add_button('text',[.02,l(3),.09,.035],'Plot levels');
add_button('edit',[.11,l(3),.06,.035],'8',[1 1 1],'levels');   
add_button('text',[.02,l(2),.09,.035],'Threshold');
add_button('edit',[.11,l(2),.06,.035],'0',[1 1 1],'t_lev');

l=.42+[0,.01,.05,.12];
a(3)=add_button('frame',[.01,l(1),.17,.04*2.5+.015]);
add_button('pushbutton',[.02,l(2),.15,.035],'Clear Markers',[],'clear',...
    [S,'''clear'');'],1,'uicontextmenu',clearMark);

add_icon([.02,l(3),.04,.05],'rect',[S,'''arect'');'],'tag','add_rect','userdata',...
    POLY,'uicontextmenu',addrect,'tooltip','Insert rectangle');
add_icon([.06,l(3),.04,.05],'poly',[S,'''poly'');'],'tag','add_poly','userdata',POLY,...
    'tooltip','Insert Polygon');
add_icon([.1,l(3),.04,.05],'text',[S,'''text'');'],'tag','add_text','userdata',P,...
    'tooltip','Add text');
add_icon([.14,l(3),.04,.05],'QPC',[S,'''qpc'');'],'tag','add_point','enable','off',...
    'tooltip','Insert a Quadratic Phase Coupling point');

l=.28+[0,.01,.05,.09,.13,.17,.21];
a(4)=add_button('frame',[.01,l(1),.17,.04*2.5+.015]);
add_icon([.02,l(3),.04,.05],'zoomin',[S,'''zoomIn'');'],'uicontextmenu',zoom_restore,...
    'tag','zoomIn');
add_icon([.07,l(3),.04,.05],'zoomout',[S,'''zoomOut'');'],'tag','zoomOut');

add_button('pushbutton',[.02,l(2),.15,.035],'Change Axes limits',[],'axesLimit',...
    [S,'''chaxes'');']);
add_button('edit',[.25,.02,.6,.05],'',[],'info','',[],'fontsize',9,...
    'horizontalalignment','left');

add_button('pushbutton',[.02,.01,.15,.06],'Draw Now',[0 1 1],'',[S,'''Refresh'');']',[],...
    'foregroundcolor',[1 1 1]);
l=.1+[0,.01,.05,.09,.13,.17,.21];
add_button('pushbutton',[.02,l(2),.15,.035],'Export',[0 1 1],'','export_fig(get(gcbo,''parent''));',[],...
    'foregroundcolor',[1 1 1]);
add_button('pushbutton',[.02,l(3),.15,.035],'Integrals',[],'',...
    [S,'''wavInt'');'],[],'uicontextmenu',wavIntegral);
[m,n]=size(FileN);
if m>0
    add_button('text',[.25,.96-(m-1)*0.03,.6,.03*m],FileN,[],'FileName');
end

h=findobj(fig2,'type','uicontrol');
h1=findobj(h,'style','edit');
set(h,'backgroundcolor',[.1529,.0824,.5412],'foregroundcolor',[1 1 1],'fontsize',9); 
set(h1,'backgroundcolor',[1,1,1],'foregroundcolor',[0 0 1]);
set(a,'backgroundcolor',[.251,.502,.502],'foregroundcolor',[.251,.502,.502]);

% Upper toolbar
add_icon([.01,.94,.04,.05],'print','printdlg(gcbf);');
add_icon([.05,.94,.04,.05],'help','helpdlg(help(''plot_3d''));');
add_icon([.09,.94,.04,.05],'save',[S,'''save'');']);
a=add_icon([.13,.94,.04,.05],'exit',[S,'''exit'');']);
set(a,'tag','close');
a=add_icon([.145,.94,.04,.05],'info',[S,'''information'');']);
set(a,'userdata',INFOTXT,'tag','INF');
%*************************************************************************
ht=0.76-0.03*(m+1)-.15;
axes('position',[.86,.15,.03,ht],'tag','hbar');
axes('position',[.25,.15,.6,ht],'tag','3dplot','drawmode','fast');
if time_series==1
    axes('position',[.25,0.76-0.03*m,.6,.2],'tag','2dplot','color',[0 0 1],'xtick',[]);
end
uimenu(umenu,'label','ln','callback',...
    'set(findobj(gcbf,''tag'',''log''),''userdata'',''loge'');');
uimenu(umenu,'label','log10','callback',...
    'set(findobj(gcbf,''tag'',''log''),''userdata'',''log10'');');
uimenu(umenu,'label','log2','callback',...
    'set(findobj(gcbf,''tag'',''log''),''userdata'',''log2'');');

uimenu(addrect,'label','Enter rectangle limits','callback',...
    [S,'''srect'');']);


p='check_uicm;';
uimenu(clearMark,'label','Clear all markers','callback',['set(findobj(gcbf,''tag'',''clear''),''userdata'',1);'...
        ,S,'''clear'');']);
uimenu(clearMark,'label','Clear rectangles','callback',[p,'set(findobj(gcbf,''tag'',''clear''),''userdata'',2);',S,'''clear'');']);
uimenu(clearMark,'label','Clear polygons','callback',[p,'set(findobj(gcbf,''tag'',''clear''),''userdata'',3);',S,'''clear'');']);
uimenu(clearMark,'label','Clear points','callback',[p,'set(findobj(gcbf,''tag'',''clear''),''userdata'',4);',S,'''clear'');']);
uimenu(clearMark,'label','Clear text','callback',[p,'set(findobj(gcbf,''tag'',''clear''),''userdata'',5);',S,'''clear'');']);

uimenu(wavIntegral,'label','Enter regions','callback',[S,'''awavInt'');']);

uimenu(zoom_restore,'label','Restore','callback',[S,'''zoom restore'');']);
uimenu(frac,'label','*F','callback',[p,'set(findobj(gcbf,''tag'',''fFactor''),''userdata'',''f'');'],'checked','on');
uimenu(frac,'label','*F^2','callback',[p,'set(findobj(gcbf,''tag'',''fFactor''),''userdata'',''f2'');']);
uimenu(frac,'label','*F^3','callback',[p,'set(findobj(gcbf,''tag'',''fFactor''),''userdata'',''f3'');']);

uimenu(filMenu,'label','3x3','callback',[p,'set(gco,''userdata'',[3,3]);']);
uimenu(filMenu,'label','5x5','callback',[p,'set(gco,''userdata'',[5,5]);']);
uimenu(filMenu,'label','7x7','callback',[p,'set(gco,''userdata'',[7,7]);'],'checked','on');
uimenu(filMenu,'label','3x1','callback',[p,'set(gco,''userdata'',[3,1]);']);
uimenu(filMenu,'label','5x1','callback',[p,'set(gco,''userdata'',[5,1]);']);
uimenu(filMenu,'label','7x1','callback',[p,'set(gco,''userdata'',[7,1]);']);
uimenu(filMenu,'label','Other','callback',[p,S,'''Filter'');']);
%************************************************************************************   
%   Select Xlim, Ylim
%************************************************************************************
function [Fx,Fy]=SelectAxesLimits(p3d);

xl=sort(get(p3d,'xlim'));
yl=sort(get(p3d,'ylim'));

xl=xl(:);yl=yl(:);
string=['a=real(get(gcbf,''currentcharacter''));if a==13,set(gcbf,''userdata'',''cont'');end;'];
fig3=figure('name','Frequency limits','units','normalized','position',[.37,.4,...
        .35,.2],'menu','none','number','off','color',[.75,.75,.75],...
    'windowstyle','modal','keypressfcn',string);
a=add_button('frame',[.05,.25,.9,.7]);
set(a,'backgroundcolor',[.251,.502,.502],'foregroundcolor',[.251,.502,.502]);
h(1)=add_button('text',[.1,.7,.3,.2],'X limits');
s='if get(gco,''userdata'')==0,set(gco,''userdata'',1);else,set(gcbf,''userdata'',''cont'');end';
M(1)=add_button('edit',[.4,.7,.2,.2],num2str(round(xl(1))),[1 1 1],'Fx',s,0);
M(2)=add_button('edit',[.7,.7,.2,.2],num2str(round(xl(2))),[1 1 1],'Fx',s,0);   
h(2)=add_button('text',[.1,.3,.3,.2],'Y limits');
M(3)=add_button('edit',[.4,.3,.2,.2],num2str(yl(1),3),[1 1 1],'Fy',s,0);
M(4)=add_button('edit',[.7,.3,.2,.2],num2str(yl(2),3),[1 1 1],'Fy',s,0);  
STR=['Fy=findobj(gcbf,''tag'',''Fy'');if get(gco,''value'')==1,S=''off'';'...
        'else;S=''on'';end,set(Fy,''enable'',S);'];
k=add_button('checkbox',[.1,.55,.8,.1],'Same Freq limits for both axes',[.4 .2 1],...
    '',STR);
if xl==yl,   
    set(k,'value',1);
    set(M(3:4),'enable','off');
else,set(k,'value',0);end

set(h,'backgroundcolor',[0.4 0.2 1],'foregroundcolor',[1 1 1],'fontsize',9); 
h1(1)=add_button('pushbutton',[.51,.01,.4,.2],'OK',[],'cont','set(gcbf,''userdata'',''cont'');');
h1(2)=add_button('pushbutton',[.09,.01,.4,.2],'Cancel',[],'','delete(gcbf);');
set(h1,'backgroundcolor',[.1529,.0824,.5412],'foregroundcolor',[1 1 1],'fontsize',9); 
set(M,'buttondownfcn','set(gco,''userdata'',0)');
waitfor(fig3,'userdata','cont');
if ~ishandle(fig3)
    Fx=[];Fy=[];return;
end
Fx=get(M(1:2),'string');
Fx=sort([str2num(Fx{1}),str2num(Fx{2})]);
if get(k,'value')==1
    Fy=Fx;
else
    Fy=get(M(3:4),'string');
    Fy=sort([str2num(Fy{1}),str2num(Fy{2})]);  
end
delete(fig3);

%************************************************************************************
%************************************************************************************   
function [x,y]=SelectRectLimits;

x=[];y=[];
fig4=figure('name','Rectangle limits','units','normalized','position',[.37,.4,...
        .3,.2],'menu','none','number','off','color',[.75,.75,.75],...
    'windowstyle','modal');
a=add_button('frame',[.05,.25,.9,.7]);
set(a,'backgroundcolor',[.251,.502,.502],'foregroundcolor',[.251,.502,.502]);
h(1)=add_button('text',[.1,.7,.3,.2],'X Interval');
M(1)=add_button('edit',[.4,.7,.2,.2],'',[1 1 1],'Fx');
M(2)=add_button('edit',[.7,.7,.2,.2],'',[1 1 1],'Fx');   
h(2)=add_button('text',[.1,.3,.3,.2],'Y Interval');
M(3)=add_button('edit',[.4,.3,.2,.2],'',[1 1 1],'Fy');
M(4)=add_button('edit',[.7,.3,.2,.2],'',[1 1 1],'Fy');  

set(h,'backgroundcolor',[0.4 0.2 1],'foregroundcolor',[1 1 1],'fontsize',9); 
h1(1)=add_button('pushbutton',[.51,.01,.48,.2],'OK',[],'cont',...
    'set(gcbf,''userdata'',''cont'');');
h1(2)=add_button('pushbutton',[.01,.01,.48,.2],'Cancel',[],'','delete(gcbf);');
set(h1,'backgroundcolor',[.1529,.0824,.5412],'foregroundcolor',[1 1 1],'fontsize',9); 
waitfor(fig4,'userdata','cont');
if ~ishandle(fig4),return,end
Fx=get(M(1:2),'string');
x=sort([str2num(Fx{1}),str2num(Fx{2})]);
Fy=get(M(3:4),'string');
y=sort([str2num(Fy{1}),str2num(Fy{2})]);  
delete(fig4);

%************************************************************************************
function x=selectZlim(minP,maxP);
x=[];
fig5=figure('name','3D plot limits','units','normalized','position',[.37,.4,...
        .3,.2],'menu','none','number','off','color',[.75,.75,.75],...
    'windowstyle','modal');
a=add_button('frame',[.05,.25,.9,.7]);
set(a,'backgroundcolor',[.251,.502,.502],'foregroundcolor',[.251,.502,.502]);
h(1)=add_button('text',[.1,.7,.5,.2],'Minimum level');
M(1)=add_button('edit',[.6,.7,.3,.2],num2str(minP,3),[1 1 1],'minP');
h(2)=add_button('text',[.1,.48,.5,.2],'Maximum level');
M(2)=add_button('edit',[.6,.48,.3,.2],num2str(maxP,3),[1 1 1],'maxP');
cb=add_button('check',[.1,.3,.8,.15],'Autoscale');

set(h,'backgroundcolor',[0.4 0.2 1],'foregroundcolor',[1 1 1],'fontsize',9); 
h1(1)=add_button('pushbutton',[.51,.01,.48,.2],'OK',[.8 .6 1],'cont',...
    'set(gcbf,''userdata'',''cont'');');
h1(2)=add_button('pushbutton',[.01,.01,.48,.2],'Cancel',[.8 .6 1],'','delete(gcbf);');
set(h1,'backgroundcolor',[.1529,.0824,.5412],'foregroundcolor',[1 1 1],'fontsize',9); 
waitfor(fig5,'userdata','cont');
if ~ishandle(fig5),return,end
if get(cb,'value')==1
    x=[];
else
    x(1)=str2num(get(M(1),'string'));
    x(2)=str2num(get(M(2),'string'));
end
delete(fig5);

%************************************************************************************
function [LF,HF]=addWavIntFig;
LF=[];HF=[];1
fig6=figure('name','Wavelets integrals','units','normalized','position',[.37,.4,...
        .26,.2],'menu','none','number','off','color',[.75,.75,.75],...
    'windowstyle','modal');
a=add_button('frame',[.05,.25,.9,.7]);
set(a,'backgroundcolor',[.251,.502,.502],'foregroundcolor',[.251,.502,.502]);

h(1)=add_button('text',[.1,.7,.3,.2],'LF region');
M(1)=add_button('edit',[.4,.7,.2,.2],'0.02',[1 1 1],'LF');
M(2)=add_button('edit',[.7,.7,.2,.2],'0.15',[1 1 1],'LF');   
h(2)=add_button('text',[.1,.3,.3,.2],'HF region');
M(3)=add_button('edit',[.4,.3,.2,.2],'0.15',[1 1 1],'HF');
M(4)=add_button('edit',[.7,.3,.2,.2],'0.5',[1 1 1],'HF');  

set(h,'backgroundcolor',[0.4 0.2 1],'foregroundcolor',[1 1 1],'fontsize',9); 
h1(1)=add_button('pushbutton',[.51,.01,.48,.2],'OK',[],'cont',...
    'set(gcbf,''userdata'',''cont'');');
h1(2)=add_button('pushbutton',[.01,.01,.48,.2],'Cancel',[],'','delete(gcbf);');
set(h1,'backgroundcolor',[.1529,.0824,.5412],'foregroundcolor',[1 1 1],'fontsize',9); 
waitfor(fig6,'userdata','cont');
if ~ishandle(fig6),return,end
x=get(M(1:2),'string');
LF=sort([str2num(x{1}),str2num(x{2})]);
y=get(M(3:4),'string');
HF=sort([str2num(y{1}),str2num(y{2})]);  
delete(fig6);

%************************************************************************************
function add_current_location(fig2,x,y,info);
if isempty(x),return;end
[P,ff,tt,series,p3d,p2d]=get_figure_parameters(fig2);

S=['plot_3d(',num2str(fig2),','];

delete(findobj(fig2,'tag','temp_point'));
[m,n]=min(abs(tt-x));
if ~isempty(series)
    [mi,ma]=min_max(series(:,1));
    line(x*[1,1],[mi,ma],'color',[0 0 0],'tag','temp_point','linewidth',1.5,'parent',p2d);    
end
set(fig2,'currentaxes',p3d);
if flag==0
    line([x,x],[ff(1),ff(end)],'color',[0 0 0],'tag','temp_point','linewidth',1.5,'parent',p3d);
    line([tt(1),tt(end)],[y,y],'color',[0 0 0],'tag','temp_point','linewidth',1.5,'parent',p3d);
    string=['Current x point: ',num2str(tt(n),5),'  ,   ',num2str(y)];
    set(info,'string',string)   
    return
end
%************************************************************************************
function add_QPC(fig2,x,y,info);
if isempty(x),return;end
[P,ff,tt,series,p3d,p2d]=get_figure_parameters(fig2);

S=['plot_3d(',num2str(fig2),','];

delete(findobj(fig2,'tag','temp_point'));
[m,n]=min(abs(tt-x));
[m1,n1]=min(abs(ff-y));
V=abs(P(n1,n));
phase=angle(P(n1,n))*180/pi;
string=['Current point (x,y): ',num2str(x,3),'  ,  ',num2str(y,3),',  Value: ',...
        num2str(V,3),'   Phase: ',num2str(phase,3)];
set(info,'string',string)   
if ~isempty(series)
    [mi,ma]=min_max(series(:,1));
    l(1)=line(x*[1,1],[mi,ma],'color',[0 0 0],'linewidth',1.5,'parent',p2d);    
    l(2)=line(y*[1,1],[mi,ma],'color',[0 0 0],'linewidth',1.5,'parent',p2d);
    l(3)=line((x+y)*[1,1],[mi,ma],'color',[0 0 0],'linewidth',1.5,'parent',p2d);
end
delete(findobj(p3d,'type','text','tag','temp_point'));
a=1.5*(tt(2)-tt(1));
aHeadXdata=x+[a,0,a];
aHeadYdata=y+[a,0,-a];
l(4)=line(aHeadXdata,aHeadYdata,'color',[1 0 0],'linewidth',2,'parent',p3d);
l(5)=line(x+[0,5*a],y+[0,0],'color',[1 0 0],'linewidth',2,'parent',p3d);
set(l,'tag','temp_point','Hittest','off')

%************************************************************************************   
function draw_line(fig2,handle,tag)
S=['plot_3d(',num2str(fig2),','];
if strmatch(get(handle,'enable'),'off'),return;end
POLY=get(handle,'userdata');
if get(findobj(fig2,'tag','drawMode'),'value')==2
    col=[0 0 0];
else,
    col=[1 0 0];
end
h=[];
for i=1:length(POLY.PEAK)
    pos=POLY.POS{i};
    peak=POLY.PEAK{i};
    MAX=POLY.MAX{i};
    Phase=POLY.Phase{i};
    VOL=POLY.VOL{i};
    rec_men=uicontextmenu;
    h(i)=line(pos.x,pos.y,'linewidth',2,'ButtonDownFcn',[S,'''inf'');'],...
        'userdata',[MAX,Phase,VOL,peak],'color',col,'tag',tag,...
        'uicontextmenu',rec_men);
  
    uimenu(rec_men,'label',['Xlim:',num2str(pos.x(1),3),' - ',num2str(pos.x(2),3)]);
    uimenu(rec_men,'label',['Ylim:',num2str(pos.y(1),3),' - ',num2str(pos.y(3),3)]);
end
%************************************************************************************   
function add_prev_text(fig2,add_text);
P=get(add_text,'userdata');
if get(findobj(fig2,'tag','drawMode'),'value')==2
    col=[0 0 0];
else,
    col=[1 0 0];
end
for k=1:length(P.TEXT)
    a=P.POS{k};te=P.TEXT{k};
    c=text(a(1),a(2),te);
    set(c,'color',col,'tag','temp_text');
end

%************************************************************************************   
function [MAX,Phase,VOL,peak]=comp_rectangle(fig2,pos);

[P,y,x,series,p3d,p2d]=get_figure_parameters(fig2);

BC=abs(P);
x1=pos.x;y1=pos.y;
x=x(:)';y=y(:)';
[mi,ma]=min_max(x1);
[mi1,ma1]=min_max(y1);
ax=find(x>mi &x<ma);ax=ax(:)';
ax=[max([1,ax(1)-1]),ax,min([length(x),ax(end)+1])];
ay=find(y>mi1 &y<ma1);ay=ay(:)';
ay=[max([1,ay(1)-1]),ay,min([length(y),ay(end)+1])];
n1=find(x>=min(ax) &x<=max(ax));
n2=find(y>=min(ay) &y<=max(ay));
f1=y(ay)';
t1=x(ax);
dx=diff(x(ax));dy=diff(y(ay));
dxy=dx'*dy;
dxy(:,end+1)=dxy(:,end);
dxy(end+1,:)=dxy(end,:);
peak=[(ma+mi)/2,(ma1+mi1)/2];
K=BC(ay,ax);K1=P(ay,ax);
MAX=max(max(K));
[m,n]=max(K(:));
Phase=angle(K1(n))*180/pi;
VOL=sum(sum(dxy'.*K));
%VOL=sum(sum(dxy'.*K))/length(K(:))
%VOL=sum(sum(dxy'.*K))/((ax(end)-ax(1))*(ay(end)-ay(1)))
VOL=sum(sum(dxy'.*K))/(ax(end)-ax(1))*60;
%VOL=sum(sum(dxy'.*K))/(ax(end)-ax(1));
%************************************************************************************   
%************************************************************************************   
function [MAX,Phase,VOL,peak]=comp_polygon(fig2,pos);
[P,y,x,series,p3d,p2d]=get_figure_parameters(fig2);

BC=abs(P);

x1=pos.x;y1=pos.y;
x=x(:)';y=y(:)';
[mi,ma]=min_max(x1);
[mi1,ma1]=min_max(y1);
ax=find(x>mi &x<ma);ax=ax(:)';
ax=[max([1,ax(1)-1]),ax,min([length(x),ax(end)+1])];
ay=find(y>mi1 &y<ma1);ay=ay(:)';
ay=[max([1,ay(1)-1]),ay,min([length(y),ay(end)+1])];
f1=y(ay)'*ones(1,length(ax));
t1=ones(length(ay),1)*x(ax);
f=inpolygon(t1,f1,x1,y1);
dx=diff(x(ax));dy=diff(y(ay));
dxy=dx'*dy;
dxy(:,end+1)=dxy(:,end);
dxy(end+1,:)=dxy(end,:);
peak=[(ma+mi)/2,(ma1+mi1)/2];
K=BC(ay,ax);K1=P(ay,ax);
MAX=max(max(K));
[m,n]=max(K(:));
Phase=angle(K1(n))*180/pi;
VOL=sum(sum(f.*dxy'.*K));

%************************************************************************************
function save_data(fig2);
[f,p]=uiputfile('*.mat','Choose output filename');
if p>0,
    [P,ff,tt,series,p3d,p2d]=get_figure_parameters(fig2);
    INFO=get(findobj(fig2,'tag','INF'),'userdata');
    figureName=get(fig2,'name');
    if isempty(series),save([p f],'P','ff','tt','figureName','INFO'),
    else ,save([p f],'P','ff','tt','series','figureName','INFO'),end  
end
%************************************************************************************
function exit_quest(fig2);

b=questdlg('Would you like to save the data?',' ','Yes','No','Cancel','Yes');
if strcmp(b,'Yes')
    [f,p]=uiputfile('*.mat','Choose output filename');
    if p>0,
        [yy,ff,tt,series,p3d,p2d]=get_figure_parameters(fig2);
        INFO=get(findobj(fig2,'tag','INF'),'userdata');
        figureName=get(fig2,'name');
        SR=1/(tt(2)-tt(1));
        if isempty(series),
            save([p f],'yy','ff','tt','figureName','INFO'),
        else ,
            save([p f],'yy','ff','tt','series','figureName','INFO','SR'),
        end  
    end
    delete(findobj(fig2,'type','axes'));close(fig2);
elseif strcmp(b,'No')
    delete(findobj(fig2,'type','axes'));close(fig2);
end 

function [P,ff,tt,series,p3d,p2d]=get_figure_parameters(fig2);
CLOSE=findobj(fig2,'tag','close');
U=get(CLOSE,'userdata');
P=U{1};ff=U{2};tt=U{3};series=U{4};
p3d=findobj(fig2,'tag','3dplot');
p2d=findobj(fig2,'tag','2dplot');