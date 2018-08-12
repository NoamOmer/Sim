function show_integrals_fcn(action,data,ax)
x=findobj(gcf,'tag','X Zoom');
y=findobj(gcf,'tag','Y Zoom');
D=35;
fig=get(x,'parent');
a=get_integrals_axes;
as=findobj(fig,'tag','Axes Signal');

l=get_integrals_handles;
switch action
case 'set Log'
    set_lines_to_log;
    comp_all_param;
    update_params;
    
case 'Select Regions'
    st=get_zoom_status;
    set_zoom_status(fig,a,l,'none');
%     set(gcbf,'pointer','fullcrosshair');
    fl=get_fix_length_status;
    clear_regions_lines
    set_all_uicontrol_status(0);
    if get_fix_length_status
        rl=get_fix_length_value;
        x1=get_region_lim_from_mouse(-rl);
        x2=get_region_lim_from_mouse(rl);
        p=sort([x1;x2]);
    else
        for k=1:3
%             waitforbuttonpress
            [p(k),y]=get_points(gcbf,1,'pointer','fullcrosshair');
            for kk=1:3
                line(p(k)*[1 1],ylim(a(kk)),'color','r','linestyle',':','tag','Region Line','parent',a(kk));
            end
            line(p(k)*[1 1],ylim(as),'color','r','linestyle',':','tag','Region Line','parent',as);
        end       
        if get(findobj(gcbf,'tag','Equal Regions'),'value')==1
            p(4)=p(3)+p(2)-p(1);
        else
            [p(4),y]=get_points(gcbf,1,'pointer','fullcrosshair');
        end 
        p=sort(p);
        for kk=1:3
            line(p(4)*[1 1],ylim(a(kk)),'color','r','linestyle',':','tag','Region Line','parent',a(kk));
        end
        line(p(4)*[1 1],ylim(as),'color','r','linestyle',':','tag','Region Line','parent',as);
    end
    set_regions(p);
    comp_all_param;
    update_params;
    set_zoom_status(fig,a,l,st);
    set_all_uicontrol_status(1);
    set(gcbf,'pointer','arrow');
    set_active('Clear Regions');
    set_active('Export Regions');
    update_time_lim(p);
    
case 'Show STD'
    set(findobj(gcbf,'tag','Show STD'),'string','Hide STD','callback','show_integrals_fcn(''Hide STD'');');
    set_std_handles('show');
case 'Hide STD'
    set(findobj(gcbf,'tag','Show STD'),'string','Show STD','callback','show_integrals_fcn(''Show STD'');');
    set_std_handles('hide');
    
case 'Export Regions'
    set(0,'ShowHiddenHandles','on')
    sf=findobj(0,'name','Save Results');
    if length(sf)==1
        p=get_time_lim_from_graphs;
        name={'Left start','Left end','Right start','Right end'};
        for k=1:length(p)
            set(findobj(sf,'tag',name{k}),'string',sprintf('%5.0f',p(k)))
        end
    end
    
case 'Clear Regions'
    clear_regions_lines;
    clear_all_params;
    update_time_lim
    set_inactive('Clear Regions');
    set_inactive('Export Regions');
    
case 'set Median'
    if get_median_flag
        set_line_to_median;
    else
        set_line_to_original;
    end
    comp_all_param;
    update_params;
    
case 'update Median'
    set_median_flag(1)
    set_line_to_median;
    comp_all_param;
    update_params;
    
case 'set X Zoom'
    set(x,'value',1);
    set(y,'value',0);
    set(fig,'windowbuttondownfcn','zoom(gcbf,''xdown'');show_integrals_fcn(''Axes X Zoom'',[],gca);');
    set(a,'buttondown','');
    set(l,'buttondown','');
    
case 'set Y Zoom'
    set(x,'value',0);
    set(y,'value',1);
    if get_log_flag
        set([a;l],'buttondown','');
        set(gcf,'WindowButtonDownFcn','zoom(gcbf,''ydown'')','WindowButtonUpFcn','ones;');
    else           
        set(a,'buttondown','show_integrals_fcn(''Y Zoom'');');
        set(l,'buttondown','show_integrals_fcn(''Y Zoom'');');
        set(fig,'windowbuttondownfcn','','WindowButtonUpFcn','');
    end
    
case 'Axes X Zoom'
    set(findobj(fig,'type','axes'),'xlim',get(ax,'xlim'));
case 'Y Zoom'
    c=get(fig,'selectiontype');
    Wh=gcbo;
    if strcmp(get(Wh,'type'),'line')
        Wh=get(Wh,'parent');
    end
    switch c
    case 'normal'
        W=get(Wh,'currentpoint');
        W=W(2,2);
        set(Wh,'ylim',[0 W]);
    case 'alt'
        ll=findobj(Wh,'type','line');
        set(Wh,'ylim',[0 max(get(ll(1),'ydata'))]);
    end
    
case 'exit'
    b=questdlg('Would you like to save the data?',' ','Yes','No','Cancel','Yes');
    if strcmp(b,'Yes')
        [f,p]=uiputfile('*.mat','Select an output file');
        if p>0,
            a=findobj(gcbf,'tag','Line Signal');
            Sig(:,1)=get(a,'xdata')';Sig(:,2)=get(a,'ydata')';
            a=findobj(gcbf,'tag','Line Int LF');
            t=get(a,'xdata');
            LF=get(a,'ydata');
            a=findobj(gcbf,'tag','Line Int HF');
            HF=get(a,'ydata');
            a=findobj(gcbf,'tag','Line Int Ratio');
            Ratio=get(a,'ydata');
            save([p f],'Sig','t','LF','HF','Ratio'),
        end
        delete(findobj(gcf,'type','axes'));close(gcf);
    elseif strcmp(b,'No')
        delete(findobj(gcf,'type','axes'));close(gcf);
    end 
end    

function g=get_median_flag
g=get(findobj(gcbf,'tag','Median Filter'),'value');

function set_median_flag(g)
set(findobj(gcbf,'tag','Median Filter'),'value',g);

function x=get_log_flag
x=get(findobj(gcbf,'tag','Log'),'value');

function D=get_median_length
D=str2num(get(findobj(gcbf,'tag','Median length'),'string'));

function l=get_integrals_handles
l(1)=findobj(gcf,'tag','Line Int LF');
l(2)=findobj(gcf,'tag','Line Int HF');
l(3)=findobj(gcf,'tag','Line Int Ratio');


function set_line_to_median
D=get_median_length;
l=get_integrals_handles;
if get_log_flag
    Q=3;
else
    Q=2;
end
for k=1:2
    y(k,:)=get(get(l(k),'parent'),'userdata');
    set(l(k),'ydata',cmed(y{k,Q},D));
end
k1=cmed(y{2,Q},D);
k2=cmed(y{1,Q},D);
k=find(abs(k1)>100*eps);
kk=find(abs(k1)<100*eps);
k3=zeros(size(k1));
k3(k)=k2(k)./k1(k);
k3(kk)=min(k3);

%set(l(3),'ydata',cmed(y{1,Q},D)./cmed(y{2,Q},D))
set(l(3),'ydata',k3)
function set_line_to_original(a,l);
l=get_integrals_handles;
if get_log_flag
    Q=3;
else
    Q=2;
end
for k=1:length(l)
    y=get(get(l(k),'parent'),'userdata');
    set(l(k),'ydata',y{Q});
end

function st=get_zoom_status;
if get(findobj(gcbf,'tag','X Zoom'),'value')
    st='x';
else
    st='y';
end

function set_zoom_status(fig,a,l,st)
switch st
case 'x'
    set(fig,'windowbuttondownfcn','zoom(gcbf,''xdown'');show_integrals_fcn(''Axes X Zoom'',[],gca);');
    set(a,'buttondown','');
    set(l,'buttondown','');
case 'y'    
    set(a,'buttondown','show_integrals_fcn(''Y Zoom'');');
    set(l,'buttondown','show_integrals_fcn(''Y Zoom'');');
    set(fig,'windowbuttondownfcn','');
case 'none'
    set(a,'buttondown','');
    set(l,'buttondown','');
    set(fig,'windowbuttondownfcn','');
end

function set_all_uicontrol_status(s);
h=findobj(gcbf,'type','uicontrol');
if s==0
    set(h,'enable','off')
else
    set(h,'enable','on')
end

function set_active(s)
set(findobj(gcbf,'tag',s),'enable','on');

function set_inactive(s)
set(findobj(gcbf,'tag',s),'enable','off');

function clear_regions_lines
h=findobj(gcbf,'tag','Region Line');
delete(h);

function set_measurement_value(s,v)
sh=findobj(gcbf,'tag',s);
if ~isempty(sh)
    set(sh,'string',num2str(v));
    if ~isempty(findstr(s,'Change'))
        set(sh,'string',sprintf('%4.1f',v));
        if v>0
            set(sh,'foregroundco','r')
        else
            set(sh,'foregroundco','b')
        end
    end
end

function v=get_measurement_value(s)
sh=findobj(gcbf,'tag',s);
if ~isempty(sh)
    v=str2num(get(sh,'string'));
else
    v=nan;
end

function comp_all_param
p=get_regions;
if ~isempty(p)
%     [V.stat_Signal_Left V.stat_Signal_STD_Left]=get_stat_from_line('Line Signal',p(1:2));
%     [V.stat_Signal_Right V.stat_Signal_STD_Right]=get_stat_from_line('Line Signal',p(3:4));
    [M,S,pVal]=get_stat_from_line('Line Signal',p);
    V.stat_Signal_Left=M(1);    V.stat_Signal_STD_Left=S(1);
    V.stat_Signal_Right=M(2);   V.stat_Signal_STD_Right=S(2);
    disp(['Signal difference: p=',num2str(pVal)])
    V.stat_Signal_Change=100*(V.stat_Signal_Right/V.stat_Signal_Left-1);
    V.stat_Signal_STD_Change=100*(V.stat_Signal_STD_Right/V.stat_Signal_STD_Left-1);
    
%     [V.stat_LF_Left V.stat_LF_STD_Left]=get_stat_from_line('Line Int LF',p(1:2));
%     [V.stat_LF_Right V.stat_LF_STD_Right]=get_stat_from_line('Line Int LF',p(3:4));
%     
    [M,S,pVal]=get_stat_from_line('Line Int LF',p);
    V.stat_LF_Left=M(1);    V.stat_LF_STD_Left=S(1);
    V.stat_LF_Right=M(2);   V.stat_LF_STD_Right=S(2);
    disp(['LF difference: p=',num2str(pVal)])
    
    V.stat_LF_Change=100*(V.stat_LF_Right/V.stat_LF_Left-1);
    V.stat_LF_STD_Change=100*(V.stat_LF_STD_Right/V.stat_LF_STD_Left-1);
    
%     [V.stat_HF_Left V.stat_HF_STD_Left]=get_stat_from_line('Line Int HF',p(1:2));
%     [V.stat_HF_Right V.stat_HF_STD_Right]=get_stat_from_line('Line Int HF',p(3:4));
    [M,S,pVal]=get_stat_from_line('Line Int HF',p);
    V.stat_HF_Left=M(1);    V.stat_HF_STD_Left=S(1);
    V.stat_HF_Right=M(2);   V.stat_HF_STD_Right=S(2);
    disp(['HF difference: p=',num2str(pVal)])
    
    V.stat_HF_Change=100*(V.stat_HF_Right/V.stat_HF_Left-1);
    V.stat_HF_STD_Change=100*(V.stat_HF_STD_Right/V.stat_HF_STD_Left-1);
    
%     [V.stat_Ratio_Left V.stat_Ratio_STD_Left]=get_stat_from_line('Line Int Ratio',p(1:2));
%     [V.stat_Ratio_Right V.stat_Ratio_STD_Right]=get_stat_from_line('Line Int Ratio',p(3:4));
       
    [M,S,pVal]=get_stat_from_line('Line Int Ratio',p);
    V.stat_Ratio_Left=M(1);    V.stat_Ratio_STD_Left=S(1);
    V.stat_Ratio_Right=M(2);   V.stat_Ratio_STD_Right=S(2);
    disp(['Ratio difference: p=',num2str(pVal)])
    
    V.stat_Ratio_Change=100*(V.stat_Ratio_Right/V.stat_Ratio_Left-1);
    V.stat_Ratio_STD_Change=100*(V.stat_Ratio_STD_Right/V.stat_Ratio_STD_Left-1);
    
    set_stats(V);
end

function [m,v,pVal]=get_stat_from_line(l,p);
[x,y]=get_line_data(l);
ind=1;
for k=1:2:4
    px=find(x>=p(k) & x<=p(k+1));
    m(ind)=mean(y(px));
    v(ind)=std(y(px));
    Y{ind}=y(px);
    ind=ind+1;
end
[f,pVal]=ttest2(Y{1},Y{2});

function [x,y]=get_line_data(s);
y=get(findobj(gcbf,'tag',s),'ydata');
x=get(findobj(gcbf,'tag',s),'xdata');


function update_params
V=get_stats;
if ~isempty(V)
    names=fieldnames(V);
    for k=1:length(names)
        q=find(names{k}=='_');
        N=names{k}((q+1):end);
        M=N;M(N=='_')=' ';
        z=getfield(V,names{k});
        set_measurement_value(M,z);
    end
end

function set_regions(p)
set_fig_data('Regions',p);

function p=get_regions
p=get_fig_data('Regions');

function set_stats(V)
set_fig_data('Stat',V);

function V=get_stats
V=get_fig_data('Stat');

function set_fig_data(S,p)
D=get(gcbf,'userdata');
D=setfield(D,S,p);
set(gcbf,'userdata',D);

function p=get_fig_data(S)
D=get(gcbf,'userdata');
p=getfield(D,S);


function clear_all_params
V=get_stats;
names=fieldnames(V);
for k=1:length(names)
    V=setfield(V,names{k},'');
end
set_stats(V);
update_params;

function set_lines_to_log
if get_median_flag
    set_line_to_median;
else
    set_line_to_original;
end
update_axes_ylim;
update_axes_zoom;
update_region_ylim;

function update_axes_zoom
show_integrals_fcn('set Y Zoom');

function update_axes_ylim
l=get_integrals_handles;
Llim=get_fig_data('Llim');
K=5;
x=get(l(1),'xdata');
%tt=find(x>(x(1)+5/Llim(1)) & x<(max(x)-5/Llim(1)));
tt=1:length(x);
if get_log_flag
    for k=1:length(l)
        y=get(l(k),'ydata');
        set(get(l(k),'parent'),'ylim',min_max(y(tt)));
    end
else
    for k=1:length(l)
        set(get(l(k),'parent'),'ylim',[0 max(get(l(k),'ydata'))]);
    end
end    
function update_region_ylim
h=findobj(gcbf,'tag','Region Line');
for k=1:length(h)
    set(h(k),'ydata',get(get(h(k),'parent'),'ylim'));
end

function a=get_integrals_axes
a(1)=findobj(gcf,'tag','Axes Int LF');
a(2)=findobj(gcf,'tag','Axes Int HF');
a(3)=findobj(gcf,'tag','Axes Int Ratio');

function set_std_handles(S);
hh=findobj(gcbf,'type','uicontrol');
r=1;
for k=1:length(hh)
    n=get(hh(k),'tag');
    if ~isempty(findstr(n,'STD')) & ~strcmp(n,'Show STD')
        h(r)=hh(k);
        r=r+1;
    end
end

switch S
case 'show'
    set(h,'visible','on')
case 'hide'
    set(h,'visible','off')
end

function update_time_lim(p)
if nargin<1
    p=get_time_lim_from_graphs;
end
ta={'Start Left','End Left','Start Right','End Right','Diff Left','Diff Right'};
for k=1:length(ta)
    if ~isempty(p)
        p(5)=p(2)-p(1);
        p(6)=p(4)-p(3);
        set(findobj(gcbf,'tag',ta{k}),'string',sprintf('%5.0f',p(k)));
    else
        set(findobj(gcbf,'tag',ta{k}),'string','');
    end
end


function p=get_time_lim_from_graphs
hp=findobj(findobj(gcbf,'tag','Axes Signal'),'tag','Region Line');
p=[];
if ~isempty(hp)
    x=get(hp,'xdata');
    x=cat(1,x{:});
    p=sort(x(:,1));
end

function fl=get_fix_length_status
fl=get(findobj(gcbf,'tag','Fix length check'),'value');

function rl=get_fix_length_value
rl=str2double(get(findobj(gcbf,'tag','Fix length edit'),'string'));


function x=get_region_lim_from_mouse(rl)
a=findobj(gcbf,'type','axes');
set(findobj(gcbf,'tag','Fix length edit'),'userdata',rl);
% waitforbuttonpress
% point=get(gca,'CurrentPoint');
[point,y]=get_points(gcbf,1,'pointer','fullcrosshair');

p(1) = point(1,1);
p(2)=p(1)+rl;
for k=1:length(a)
    axes(a(k));
    l1(k)=line([p(1) p(1)],ylim','color','r','linestyle',':','tag','R1');
    l2(k)=line([p(2) p(2)],ylim','color','r','linestyle',':','tag','R2');
end
set(gcbf,'WindowButtonMotionFcn','p=get(gca,''CurrentPoint'');set(findobj(gcbf,''tag'',''R1''),''xdata'',[p(1,1) p(1,1)]);set(findobj(gcbf,''tag'',''R2''),''xdata'',[p(1,1) p(1,1)]+get(findobj(gcbf,''tag'',''Fix length edit''),''userdata''));');
set(gcbf,'WindowButtonUpFcn','set(gcbf,''WindowStyle'',''normal'');set(gcbf,''WindowButtonMotionFcn'','''');set(findobj(gcbf,''tag'',''Fix length edit''),''userdata'',0);');
set(gcbf,'WindowStyle','modal');
waitfor(findobj(gcbf,'tag','Fix length edit'),'userdata',0)
L1=findobj(gca,'tag','R1');
L2=findobj(gca,'tag','R2');
Q1=get(L1,'xdata');
Q2=get(L2,'xdata');
x=[Q1(1);Q2(1)];
set(findobj(gcbf,'tag','R1'),'tag','Region Line');
set(findobj(gcbf,'tag','R2'),'tag','Region Line');


    