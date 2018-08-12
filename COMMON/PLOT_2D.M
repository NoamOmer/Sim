function plot_2d(varargin)
F=[];P=[];INFOTXT='';FileN=[];
K=varargin;
if length(K)>0
    if length(K{1}(:))==1
        fig2=K{1};
        choice=K{2};
    else
        fig2=[];
    end
end

if isempty(fig2)
    l=[];
    for k=1:length(K)
        if isstr(K{k})
            if findstr(lower(K{k}),'info')
                INFOTXT=K{k+1};
                l=k;
            end
        end
    end
    if ~isempty(l),    K(l:l+1)=[];      end
    for k=1:length(K)
        if isstr(K{k})
            if findstr(lower(K{k}),'filename')
                FileN=K{k+1};
                l=k;
            end
        end
    end
    if ~isempty(l),    K(l:l+1)=[];      end
    P=varargin{1};
    K(1)=[];
    if ~isempty(K)
        if ~isstr(K{1}),F=K{1};K(1)=[];end
    end
    if isempty(K)
        figureName='plot 2D';
    else
        figureName=K{1};
    end
    fig2=add_figure(F,P,FileN,INFOTXT,figureName);
    choice=0;
else
    fig2=K{1};
    choice=K{2};
end

col=[0 0 1;1 0 0;0 1 0];
info=findobj(fig2,'tag','info');
cl=findobj(fig2,'tag','close');
X=get(cl,'userdata');
P=X{1};F=X{2};
rangeF=range(F);
maxf=max(F);
BAND=findobj(fig2,'tag','band');
band=str2num(get(BAND,'string'));
Sli=findobj(fig2,'tag','lim');
set(Sli,'enable','on');
axe=findobj(fig2,'type','axes');
if ~band>0
    band=rangeF/5;
    set(BAND,'string',num2str(band))   
elseif band>=rangeF
    band=rangeF;
    set(BAND,'string',num2str(band))   
    set(Sli,'value',F(1),'enable','off');
    for k=1:length(axe)
        YY=min_max(P(:,k));
        if range(YY)==0
            YY=YY(1)+[-0.05,0.05];
        end
        a=findobj(axe,'userdata',k);
        set(a,'xlim',[F(1),F(end)],'ylim',YY),
    end
end

switch choice
case 1
    B=get(Sli,'value');
    if B>maxf-band,
        B=maxf-band;
    end,
    if strcmp(lower(get(Sli,'enable')),'on')
        set(Sli,'value',B,'max',maxf-band,'sliderstep',[.1,.9]*band/rangeF),
    end
    l=find(F>B & F<B+band);
    yy=min_max(P(l,:));        
    for k=1:length(axe)
        YY=yy(:,k);
        if range(YY)==0
            YY=YY(1)+[-0.05,0.05];
        end
        a=findobj(axe,'userdata',k);       
        set(a,'xlim',B+[0,band],'ylim',YY),
    end
case 2 %add line
    [x,y]=get_points(gca,2);
    px=find(F>=x(1) & F<=x(2));
    pos.x=[F(px(1));F(px);F(px(end));F(px(1))];
    for k=1:length(axe)
        pos.y=[min(P(:,k));P(px,k);min(P(:,k));min(P(:,k))];
        [m,n]=max(pos.y);
        a=findobj(axe,'userdata',k);       
        u=uicontextmenu;
        patch(pos.x,pos.y,col(1,:),'userdata',[m,sum(P(px,k)),pos.x(n)],...
            'buttondownfcn',['plot_2d(',num2str(fig2),',4);'],'tag','peak','parent',a,...
            'uicontextmenu',u);
        uimenu(u,'label','Delete','callback','delete(gco)');     
        uimenu(u,'label','Change pattern','callback','changepatchpattern(gco);');     
    end
case 3 %clear lines
    delete(findobj(fig2,'tag','peak'));
    set(info,'string','');
    return
case 4 
    A=get(gco,'userdata');
    [mi,ma]=min_max(get(gco,'xdata'));
    stri=['Interval: (',num2str(mi,3),',',num2str(ma,3),'),  Max value (x,y): ',...
            num2str(A(3),3),',',num2str(A(1),3),...
            '  Integral: ',num2str(A(2),5)];
    set(info,'string',stri);
case 5
    b=questdlg('Would you like to save the spectrum data?',' ','Save as MAT','Save as TXT','No','Cancel');
    if strcmp(b,'Save as MAT')
        [f,p]=uiputfile('*.mat','Select output filename');
        if p>0
            save([p,f],'P','F')
        end      
        delete(fig2);
    elseif strcmp(b,'Save as TXT')
        [f,p]=uiputfile('*.txt','Select output filename');
        if p>0
            P=[F(:),P(:)];
            save([p,f,'.txt'],'P','-ascii','-tabs')
        end      
        delete(fig2);
    elseif strcmp(b,'No')
        delete(fig2);
    end
case 6
    x=get(gca,'currentpoint');x=x(1);
    X=get(findobj(gca,'tag','psd_data'),'xdata');
    Y=get(findobj(gca,'tag','psd_data'),'ydata');
    [k,I]=min(abs(X-x));
    stri=['(x,y): ',num2str(X(I),5),',',num2str(Y(I),5)];
    set(info,'string',stri);   
case 7
    h=msgbox(get(gcbo,'userdata'),'createMode','modal'); 
case 8 %Select line
    [LF,HF,TP]=addIntFig;
    X=[{LF},{HF},{TP}];
    for mm=1:3
        x=X{4-mm};
        if ~isempty(x)
            px=find(F>x(1) & F<x(2));
            pos.x=[F(px(1));F(px);F(px(end))];
            for k=1:length(axe)
                pos.y=[min(P(:,k));P(px,k);min(P(:,k))];
                [m,n]=max(pos.y);
                a=findobj(axe,'userdata',k);       
                patch(pos.x,pos.y,col(mm,:),'userdata',[m,sum(P(px,k))*(F(2)-F(1)),pos.x(n)],...
                    'buttondownfcn',['plot_2d(',num2str(fig2),',4);'],'tag','peak','parent',a);
            end
        end
    end    
case 9 %export
    a=findobj(gca,'type','patch');
    if isempty(a),return,end
    for k=1:length(a)
        A=get(a(k),'userdata');
        [mi,ma]=min_max(get(a(k),'xdata'));
        Interval(k,:)=[mi,ma];        
        Integral(k)=A(2);
    end    
    assignin('base','Interval',Interval)
    assignin('base','Integral',Integral(:))
end

%*************************************************************************************
function fig2=add_figure(F,P,FN,INFOTXT,figureName);
fig2=figure('tag','2D plot','menu','none','number','off','units','normalized',...
    'position',[.04,.06,.9,.9],'color',[.1,.1,.1]+.6,'name',figureName);
S=['plot_2d(',num2str(fig2),','];
col=[1 1 1];
colo='bgrykmc';
[m,n]=size(P);
if n>m,P=P';end
[m,K]=size(P);
if isempty(F);F=1:m;F=F(:);end

[l,m]=size(FN);
if l>0
    a=add_button('text',[0.2,.95-(l-1)*0.03,0.78,0.03*l],FN,[],'siginf');
    set(a,'foregroundcolor',[1 1 1],'backgroundcolor',[0.3,0.3,0.7])
end
pp=0.2*K;
dp=(0.65)/K;
ppp=0.2+0.05+(0:K)*dp;ppp=ppp(end:-1:1);ppp(1)=[];
for k=1:K
    a=axes('units','normalized','position',[.2,ppp(k),.78,dp-.05],'buttondownfcn',[S,'6);'],...
        'color',col,'userdata',k);
    line(F,P(:,k),'color',colo(k),'tag','psd_data','buttondownfcn',[S,'6);'],'parent',a);
    [m,n]=min_max(F);
    l=find(F<=m+(n-m)/5);
    yy=min_max(P(l,k));
    if range(yy)==0,yy=min(yy)+[-0.5,0.5];end
    set(a,'xlim',[F(l(1)),F(l(end))],'ylim',yy)
    if k<K;
        set(a,'xtick',[])
    end
end
pos=get(a,'position');
pos=[pos(1),.14,pos(3),.04];
an=uicontrol('Style','slider','units','normalized','min',m,'max',m+(n-m)*.95,'position',...
    pos,'sliderstep',[.1 .9]/5,'tag','lim','callback',[S,'1)'],'value',m);
add_button('frame',[.44,.07,.22,.06],'',[.251,.502,.502]);

h(1)=add_button('text',[.45 .08 .05 .04],'<---');
add_button('edit',[.5 .08 .1 .04],num2str(range(F(l))),[1 1 1],'band',[S,'1)']);
h(2)=add_button('text',[.6 .08 .05 .04],'--->');
set(h,'fontsize',10);

add_button('frame',[.01,.07,.12,.29],'',[.251,.502,.502]);
h(3)=add_button('pushbutton',[.02,.29,.1,.06],'Select region',[],'',[S,'2);']);
h(4)=add_button('pushbutton',[.02,.22,.1,.06],'Enter regions',[],'',[S,'8);']);
h(5)=add_button('pushbutton',[.02,.15,.1,.06],'Clear markers',[],'',[S,'3);']);   
h(6)=add_button('pushbutton',[.02,.08,.1,.06],'Export',[],'Export',[S,'9);']);
h(7)=add_button('push',[.01,.85,.12,.06],'Export figure',[],'','export_fig(gcf);');

set(h,'backgroundcolor',[.1529,.0824,.5412],'foregroundcolor',[1 1 1]);

a=add_button('edit',[.25,.01,.6,.04],'',[1 1 1],'info');
set(a,'horizontalalignment','left','fontsize',10,'foregroundcolor',[1 0 0]);

add_icon([.01,.94,.04,.05],'print','printdlg(gcbf);');
add_icon([.055,.94,.04,.05],'help','helpdlg(help(''plot_2d''));');
a=add_icon([.1,.94,.04,.05],'exit',[S,'5);']);
set(a,'tag','close','userdata',[{P},{F}]);
a=add_icon([.145,.94,.04,.05],'info',[S,'7);']);
set(a,'userdata',INFOTXT,'tag','INF');

function [LF,HF,TP]=addIntFig;
LF=[];HF=[];TP=[];
fig6=figure('name','PSD integrals','units','normalized','position',[.37,.4,...
        .36,.2],'menu','none','number','off','color',[.75,.75,.75],...
    'windowstyle','modal');
a=add_button('frame',[.05,.25,.9,.75]);
set(a,'backgroundcolor',[.251,.502,.502],'foregroundcolor',[.251,.502,.502]);

cb=['if get(gcbo,''value'')==0,'...
        'set(findobj(gcbf,''tag'',get(gcbo,''userdata'')),''enable'',''off''),else,'...
        'set(findobj(gcbf,''tag'',get(gcbo,''userdata'')),''enable'',''on''),end,'];
h(1)=add_button('check',[.1,.78,.3,.18],'LF region',[],'',cb,'LF','value',1);
M(1)=add_button('edit',[.4,.78,.2,.18],'0.04',[1 1 1],'LF');
M(2)=add_button('edit',[.7,.78,.2,.18],'0.18',[1 1 1],'LF');   
h(2)=add_button('check',[.1,.53,.3,.18],'HF region',[],'',cb,'HF','value',1);
M(3)=add_button('edit',[.4,.53,.2,.18],'0.2',[1 1 1],'HF');
M(4)=add_button('edit',[.7,.53,.2,.18],'0.4',[1 1 1],'HF');  
h(3)=add_button('check',[.1,.28,.3,.18],'Total region',[],'',cb,'TP','value',1);
M(5)=add_button('edit',[.4,.28,.2,.18],'0.04',[1 1 1],'TP');
M(6)=add_button('edit',[.7,.28,.2,.18],'0.6',[1 1 1],'TP');  

set(h,'backgroundcolor',[0.4 0.2 1],'foregroundcolor',[1 1 1],'fontsize',9); 
h1(1)=add_button('pushbutton',[.51,.01,.48,.2],'OK',[],'cont',...
    'set(gcbf,''userdata'',''cont'');');
h1(2)=add_button('pushbutton',[.01,.01,.48,.2],'Cancel',[],'','delete(gcbf);');
set(h1,'backgroundcolor',[.1529,.0824,.5412],'foregroundcolor',[1 1 1],'fontsize',9); 
waitfor(fig6,'userdata','cont');
if ~ishandle(fig6),return,end
if get(h(1),'value')==1
    x=get(M(1:2),'string');
    LF=sort([str2num(x{1}),str2num(x{2})]);
end
if get(h(2),'value')==1
    y=get(M(3:4),'string');
    HF=sort([str2num(y{1}),str2num(y{2})]);  
end
if get(h(3),'value')==1
    y=get(M(5:6),'string');
    TP=sort([str2num(y{1}),str2num(y{2})]);  
end
delete(fig6);
