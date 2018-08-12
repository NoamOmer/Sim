function setting_file=createNewSet(filen);
CM={'HSV','GRAY','HOT','COOL','BONE','COPPER','PINK','LAG','PRISM','JET'};
COLOR={'None','White','Blue','Yelow','Red','Magenta','Cyan','Black'};
l=0.2+[0,0.01:0.08:0.7];l=l(9:-1:1);
dy=0.06;
x=.05;
dx=.49;
%fig=figure('windowstyle','modal');
fig=figure('units','normalized','position',[0.2,0.2,0.4,0.6]);

add_button('frame',[x-0.01,l(end)-.01,1.02-2*x,l(1)-l(9)+.02+dy],'',[.251,.502,.502]);

a(1)=add_button('text',[x,l(1),dx,dy],'Paper position:');
a(2)=add_button('text',[x,l(2),dx,dy],'Colormap:');
a(3)=add_button('text',[x,l(3),dx,dy],'Axes line color:');
a(4)=add_button('text',[x,l(4),dx,dy],'Axes background color:');
a(5)=add_button('text',[x,l(5),dx,dy],'Graph line color:');
a(6)=add_button('text',[x,l(6),dx,dy],'Text color:');
a(7)=add_button('text',[x,l(7),dx,dy],'Text font size:');
a(8)=add_button('text',[x,l(8),dx,dy],'Axes font (labels) size:');

x=.55;
dx=.4;
b(1)=add_button('edit',[x,l(1),dx,dy],'[0,0,8.6,11]',[1 1 1],'PaperPos');
b(2)=add_button('check',[x,l(2),.12,dy],'1-',[],'1-map','',[],'value',0);
c(1)=add_button('popup',[x+.122,l(2),dx-.122,dy],CM);
c(2)=add_button('popup',[x,l(3),dx,dy],COLOR,[],'','',[],'value',7);
c(3)=add_button('popup',[x,l(4),dx,dy],COLOR);
c(4)=add_button('popup',[x,l(5),dx,dy],COLOR,[],'','',[],'value',7);
c(5)=add_button('popup',[x,l(6),dx,dy],COLOR);

b(3)=add_button('edit',[x,l(7),dx,dy],'10',[1 1 1],'TextSize');
b(4)=add_button('edit',[x,l(8),dx,dy],'9',[1 1 1],'AxesFontSize');


but(1)=add_button('push',[.2,.01,.25,.06],'OK',[],'','set(gcbf,''userdata'',1);');
but(2)=add_button('push',[.55,.01,.25,.06],'Cancel',[],'','set(gcbf,''userdata'',0);');

set(c,'backgroundcolor',[1,1,1])
set(a,'backgroundcolor',[.25,.2,.9],'foregroundcolor',[1 1 1],'fontsize',9); 
set(but,'backgroundcolor',[.2,.1,.7],'foregroundcolor',[1,1,1],'fontweight','bold',...
   'fontsize',12);
if nargin==1
    load(filen,'-mat')
    set(b(3),'string',num2str(S.textfontsize));
    set(b(4),'string',num2str(S.axesfontsize));
    set(b(1),'string',num2str(S.paperpos));
    set(c(2),'value',find(strcmp(lower(COLOR),lower(S.plot2dlinecolor))))
    set(c(3),'value',find(strcmp(lower(COLOR),lower(S.plot2dbackcolor))))
    set(c(4),'value',find(strcmp(lower(COLOR),lower(S.plot2dlinecolor))))   
    set(c(3),'value',find(strcmp(lower(COLOR),lower(S.textcolor))))   
else
    filen='*.set';
end

waitfor(fig,'userdata');
if isempty(fig),return,end
if get(fig,'userdata')==1
    S.paperpos=str2num(get(b(1),'string'));
    m=CM{get(c(1),'value')};
    S.map=eval(m);
    if get(b(2),'value')==1
        S.map=1-S.map;
    end
    C=COLOR{get(c(2),'value')};
    S.axeslinecolor=C;
    C=COLOR{get(c(3),'value')};
    S.plot2dbackcolor=C;
    C=COLOR{get(c(4),'value')};
    S.plot2dlinecolor=C(1);
    C=COLOR{get(c(5),'value')};
    S.textcolor=lower(C);
    S.textfontsize=str2num(get(b(3),'string'));
    S.axesfontsize=str2num(get(b(4),'string'));    
    a=questdlg('Would you like to save the new set ?','','Yes','No','Yes');
    if strcmp(a,'Yes')
        [f,p]=uiputfile(filen,'Select output file name for the new set');
        if p>1
            fn=[p,f];
            setting_file=[get_file_prop(fn,'path'),get_file_prop(fn,'name'),'.set'];
            save (setting_file,'-mat','S');
        end
    else
        setting_file=S;
    end
end
delete(fig);