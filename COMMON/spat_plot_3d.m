function spat_plot_3d(action);
if nargin==0
    action='newPlot';
end
switch action
case 'newPlot'
    [f,p]=uigetfile('*.mat','Select file to plot');
    if p==0;return;end
    S=load([p,f]);
    NA=fieldnames(S);
    INFO=[];
    tt=[];
    ff=[];
    time_series=[];
    yy=[];
    for k=1:length(NA)
        a=S.(NA{k});
        if iscell(a)
            INFO=NA{k};
        else
            [m,n]=size(a);
            if m*n>1
                if m==1 | n==1
                    if isempty(ff)
                        ff=NA{k};
                    else
                        tt=NA{k};
                    end
                elseif m==2 | n==2
                    time_series=NA{k};
                else
                    yy=NA{k};
                end
            end
        end
    end
    if ~isempty(yy) & ~isempty(tt) & ~isempty(ff)
        [m,n]=size(S.(yy));
        [m1,n1]=size(S.(ff));
        if m1*n1==m
            a=ff;
            ff=tt;
            tt=a;
        end
    end
    dx=0.08;  
    l=linspace(0.01,6*dx,7);
    l=.15+[0,l(end:-1:1)]; 
    fig=figure('units','normalized','position',[.2,.2,.3,.5]);
    a=add_button('frame',[.01,l(1),98,dx*7+0.01]);
    set(a,'backgroundcolor',[.251,.502,.502],'foregroundcolor',[.251,.502,.502]);
    dx=dx-.01;
    U.FN=add_button('text',[.01,.88,.98,dx],[p,f]);
    set(U.FN,'tag','FN','backgroundcolor',[.1529,.0824,.5412],'foregroundcolor',[1 1 1],'fontsize',9); 


    U.fields=add_button('popup',[.1,.75,.8,dx],NA,[1 1 1],'','',S);
    
    
    U.matrix=add_button('edit',[.51,l(2),.45,dx],yy,[1 1 1],'matrix','',[]);
    h(1)=add_button('pushbutton',[.04,l(2),.45,dx],'Matrix',[],'','spat_plot_3d(''Select'');',U.matrix);
    
    U.x=add_button('edit',[.51,l(3),.45,dx],tt,[1 1 1],'xvector','',[]);
    h(2)=add_button('pushbutton',[.04,l(3),.45,dx],'X vector',[],'','spat_plot_3d(''Select'');',U.x);
    
    U.y=add_button('edit',[.51,l(4),.45,dx],ff,[1 1 1],'yvector','',[]);
    h(3)=add_button('pushbutton',[.04,l(4),.45,dx],'Y vector',[],'','spat_plot_3d(''Select'');',U.y);
    
    U.time_series=add_button('edit',[.51,l(5),.45,dx],time_series,[1 1 1],'series','',[]);
    h(4)=add_button('pushbutton',[.04,l(5),.45,dx],'Time series',[],'','spat_plot_3d(''Select'');',U.time_series);

    U.info=add_button('edit',[.51,l(6),.45,dx],INFO,[1 1 1],'info','',[]);
    h(5)=add_button('pushbutton',[.04,l(6),.45,dx],'Info text',[],'','spat_plot_3d(''Select'');',U.info);
    
    U.first=add_button('edit',[.3,l(7),.2,dx],'',[1 1 1],'first','',[]);
    h(6)=add_button('text',[.04,l(7),.26,dx],'First pts',[],'','',U.first);
    
    U.last=add_button('edit',[.77,l(7),.2,dx],'',[1 1 1],'last','',[]);
    h(7)=add_button('text',[.51,l(7),.26,dx],'Last pts',[],'','spat_plot_3d(''Select'');',U.last);
    
    U.all_pts=add_button('check',[.04,l(8),.45,dx],'Use all points');
    set(U.all_pts,'value',1);
    set(fig,'userdata',U);
    
  
    g(1)=add_button('pushbutton',[.01,.01,.3,dx],'PLOT',[],'','spat_plot_3d(''Plot'');');
    g(2)=add_button('pushbutton',[.33,.01,.3,dx],'Refresh',[],'','spat_plot_3d(''Refresh'');');
    g(3)=add_button('pushbutton',[.66,.01,.3,dx],'Close',[],'','spat_plot_3d(''Close'');');
    set(h,'backgroundcolor',[.1529,.0824,.5412],'foregroundcolor',[1 1 1],'fontsize',9); 
    set(g,'backgroundcolor',[.1529,.0824,.412],'foregroundcolor',[1 1 1],'fontsize',9); 
case 'Plot'
    U=get(gcbf,'userdata');
    S=get(U.fields,'userdata');
    m=get(U.matrix,'string');
    if isempty(m),return;end
    yy=getfield(S,m);
    [mm,nn]=size(yy);
    m=get(U.x,'string');
    if isempty(m),
        xvector=1:mm;
    else
        xvector=getfield(S,m);
    end
    m=get(U.y,'string');
    if isempty(m),
        yvector=1:nn;
    else
        yvector=getfield(S,m);
    end
    m=get(U.time_series,'string');
    if isempty(m),
        s=[];
    else
        s=getfield(S,m);
    end
    figname='Plot 3d';
    m=get(U.info,'userdata');
    if isempty(m),
        INFO={''};
    else
        INFO=getfield(S,m);
    end
    k1=1;
    k2=mm;
    if get(U.all_pts,'value')==0
        k1=str2num(get(U.first,'string'));
        k2=str2num(get(U.last,'string'));
    end
    yy=yy(k1:k2,:);
    xvector=xvector(k1:k2);
    if isempty(s)
        plot_3d(yy',yvector,xvector,figname,'info',INFO);
    else
        plot_3d(yy',yvector,xvector,s,figname,'info',INFO);
    end
case 'Refresh'
    [f,p]=uigetfile('*.mat','Select file to plot');
    if p==0;return;end
    U=get(gcbf,'userdata');
    set(U.FN,'string',[p,f])
    S=load([p,f]);
    NA=fieldnames(S);
    set(U.fields,'string',NA,'userdata',S)
case 'Close'
    close(gcf)
case 'Select'
    l=findobj(gcbf,'style','popup');
    S=get(l,'string');
    k=get(l,'value');
    D=S{k};
    set(get(gcbo,'userdata'),'string',D);
end

