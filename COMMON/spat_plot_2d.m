function spat_plot_2d(action);
if nargin==0
    action='newPlot';
end
switch action
case 'newPlot'
    [f,p]=uigetfile('*.mat','Select file to plot');
    if p==0;return;end
    S=load([p,f]);
    NA=fieldnames(S);
    
    fig=figure('units','normalized','position',[.2,.2,.3,.5]);
    l=.15+[0,.6,.48,.36,.24,.12,.01];
    a=add_button('popup',[.1,.9,.8,.08],NA,[1 1 1],'','',S);
    
    a=add_button('frame',[.01,l(1),98,.12*6]);
    set(a,'backgroundcolor',[.251,.502,.502],'foregroundcolor',[.251,.502,.502]);
    
    a=add_button('edit',[.51,l(3),.45,.1],'',[1 1 1],'xvector','',[]);
    h(1)=add_button('pushbutton',[.04,l(3),.45,.1],'X vector',[],'','spat_plot_2d(''Select'');',a);
    
    a=add_button('edit',[.51,l(4),.45,.1],'',[1 1 1],'yvector','',[]);
    h(2)=add_button('pushbutton',[.04,l(4),.45,.1],'Y vector',[],'','spat_plot_2d(''Select'');',a);
    
    g(1)=add_button('pushbutton',[.01,.01,.3,.1],'PLOT',[],'','spat_plot_2d(''Plot'');');
    g(2)=add_button('pushbutton',[.33,.01,.3,.1],'Refresh',[],'','spat_plot_2d(''Refresh'');');
    g(3)=add_button('pushbutton',[.66,.01,.3,.1],'Close',[],'','spat_plot_2d(''Close'');');
    set(h,'backgroundcolor',[.1529,.0824,.5412],'foregroundcolor',[1 1 1],'fontsize',9); 
    set(g,'backgroundcolor',[.1529,.0824,.412],'foregroundcolor',[1 1 1],'fontsize',9); 
case 'Plot'
    l=findobj(gcbf,'style','popup');
    S=get(l,'userdata');
    m=get(findobj(gcbf,'tag','yvector'),'userdata');
    if isempty(m),return;end
    yvector=getfield(S,m);
    [mm,nn]=size(yvector);
    m=get(findobj(gcbf,'tag','xvector'),'userdata');
    if isempty(m),
        xvector=1:mm;
    else
        xvector=getfield(S,m);
    end
    plot_2d(yvector,xvector);
case 'Refresh'
    [f,p]=uigetfile('*.mat','Select file to plot');
    if p==0;return;end
    l=findobj(gcbf,'style','popup');
    S=load([p,f]);
    NA=fieldnames(S);
    set(l,'string',NA,'userdata',S)
case 'Close'
    close(gcf)
case 'Select'
    l=findobj(gcbf,'style','popup');
    S=get(l,'string');
    k=get(l,'value');
    D=S{k};
    set(get(gcbo,'userdata'),'string',D,'userdata',D);
end

