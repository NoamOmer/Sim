function varargout=show_integrals(tx,x,tM,f,P,Llim,Hlim,info,FileN)
f=f(:)';
[m,n]=size(P);
len=length(f);
if len==n,P=P';end
le=size(P,2);
fl=find(f>Llim(1) & f<Llim(2));
fh=find(f>Hlim(1) & f<Hlim(2));
L=(P(fl(1:(end-1)),:)+P(fl(2:end),:))'/2*diff(f(fl)');
H=(P(fh(1:(end-1)),:)+P(fh(2:end),:))'/2*diff(f(fh)');
Total=(P(1:(end-1),:)+P(2:end,:))'/2*diff(f');
info=[];
if ~isempty(info)
    ti{1}=info{1};
    ti{2}=[info{2} ' , ' info{5} ' , ' info{6}];
else
    ti='The signal';
end

if nargout==0
    NS=4;
    d1=(0.95-.11)/NS;
    d2=0.0527;
    d3=d1-d2;
    d4=.4;
    d5=.95-d4;
    p=.11+((NS-1):-1:0)*d1;
    fig=figure('menu','none','unit','norm','pos',[0.02    0.05    0.97    0.91],'name','Integrals','Doublebuffer','on');
    S.type='Integrals';
    S.Regions=[];
    S.Stat=[];
    S.Hlim=Hlim;
    S.Llim=Llim;
    set(fig,'userdata',S);
    orient tall;
    for k=1:NS
        ax(k)=axes('unit','norm','pos',[d4 p(k) d5 d3],'drawmode','fast');
    end
    add_button('text',[d4,.94,d5,.04],FileN);
    LH=L./H;
    axes(ax(1))
    set(gca,'userdata',{tx,x})
    set(gca,'tag','Axes Signal');
    line(tx,x,'tag','Line Signal')
    title(ti,'interpre','none')
    axis tight
    set(gca,'xticklabel','');
    m1=5*floor(min(x)/5);
    m2=5*ceil(max(x)/5);
    % ylim([m1 m2]);
    axes(ax(2))
    set(gca,'userdata',{tM,L,10*log10(L)})
    set(gca,'tag','Axes Int LF');
    line(tM,L,'tag','Line Int LF');
    axis tight
    title(sprintf('LF peak [%2.2f-%2.2f Hz]',Llim));
    set(gca,'xticklabel','');
    y=ylim;
    %  ylim([0 y(2)]);
    axes(ax(3))
    set(gca,'userdata',{tM,H,10*log10(H)})
    set(gca,'tag','Axes Int HF');
    line(tM,H,'tag','Line Int HF');
    axis tight
    title(sprintf('HF peak [%2.2f-%2.2f Hz]',Hlim));
    set(gca,'xticklabel','');
    y=ylim;
    %    ylim([0 y(2)]);
    axes(ax(4))
    set(gca,'userdata',{tM,LH,10*log10(LH)})
    set(gca,'tag','Axes Int Ratio');
    line(tM,LH,'tag','Line Int Ratio');
    title('LF/HF');
    xlabel('Time')
    y=ylim;
    axis tight
    % ylim([0 y(2)]);
    %    set(a,'buttondownfcn',['set(findobj(gcf,''type'',''axes''),''ylimmode'',''auto'',''xlim'',get(gca,''xlim''));']);
    hig=0.03;
    k1=1;
    ho1=0.04;
    ho2=ho1+0.005;
    l=0.03;
    q1=0.02;
    q2=d4-0.08;
    q3=.07;
    offs=0.005;
    q4=q1+q3+offs;
    q5=(q2+q1-q4-2*offs)/3;
    q6=q4+q5+offs;
    q7=q6+q5+offs;
    
    h(k1)=add_button('text',[q1,l,q3,hig],'Ratio std:',[0.4 0.2 1],'Ratio STD text','','','foregroundcolor',[1 1 1],'fontsize',9,'visible','off');
    k1=k1+1;
    h(k1)=add_button('edit',[q4,l,q5,hig],'',[0.4 0.2 1],'Ratio STD Left','','','HorizontalAlignment','left','foregroundcolor',[1 1 1],'fontsize',9,'visible','off');
    k1=k1+1;
    h(k1)=add_button('edit',[q6,l,q5,hig],'',[0.4 0.2 1],'Ratio STD Right','','','HorizontalAlignment','left','foregroundcolor',[1 1 1],'fontsize',9,'visible','off');
    k1=k1+1;
    h(k1)=add_button('edit',[q7,l,q5,hig],'','w','Ratio STD Change','','','HorizontalAlignment','left','foregroundcolor',[1 1 1],'fontsize',9,'visible','off');
    k1=k1+1;l=l+ho2;
    h(k1)=add_button('text',[q1,l,q3,hig],'Ratio mean:',[0.4 0.2 1],'','','','foregroundcolor',[1 1 1],'fontsize',9);
    k1=k1+1;
    h(k1)=add_button('edit',[q4,l,q5,hig],'',[0.4 0.2 1],'Ratio Left','','','HorizontalAlignment','left','foregroundcolor',[1 1 1],'fontsize',9);
    k1=k1+1;
    h(k1)=add_button('edit',[q6,l,q5,hig],'',[0.4 0.2 1],'Ratio Right','','','HorizontalAlignment','left','foregroundcolor',[1 1 1],'fontsize',9);
    k1=k1+1;
    h(k1)=add_button('edit',[q7,l,q5,hig],'','w','Ratio Change','','','HorizontalAlignment','left','foregroundcolor',[1 1 1],'fontsize',9);
    k1=k1+1;l=l+ho2;
    h(k1)=add_button('text',[q1,l,q3,hig],'HF std:',[0.4 0.2 1],'HF STD text','','','foregroundcolor',[1 1 1],'fontsize',9,'visible','off');
    k1=k1+1;
    h(k1)=add_button('edit',[q4,l,q5,hig],'',[0.4 0.2 1],'HF STD Left','','','HorizontalAlignment','left','foregroundcolor',[1 1 1],'fontsize',9,'visible','off');
    k1=k1+1;
    h(k1)=add_button('edit',[q6,l,q5,hig],'',[0.4 0.2 1],'HF STD Right','','','HorizontalAlignment','left','foregroundcolor',[1 1 1],'fontsize',9,'visible','off');
    k1=k1+1;
    h(k1)=add_button('edit',[q7,l,q5,hig],'','w','HF STD Change','','','HorizontalAlignment','left','foregroundcolor',[1 1 1],'fontsize',9,'visible','off');
    k1=k1+1;l=l+ho2;
    h(k1)=add_button('text',[q1,l,q3,hig],'HF mean:',[0.4 0.2 1],'','','','foregroundcolor',[1 1 1],'fontsize',9);
    k1=k1+1;
    h(k1)=add_button('edit',[q4,l,q5,hig],'',[0.4 0.2 1],'HF Left','','','HorizontalAlignment','left','foregroundcolor',[1 1 1],'fontsize',9);
    k1=k1+1;
    h(k1)=add_button('edit',[q6,l,q5,hig],'',[0.4 0.2 1],'HF Right','','','HorizontalAlignment','left','foregroundcolor',[1 1 1],'fontsize',9);
    k1=k1+1;
    h(k1)=add_button('edit',[q7,l,q5,hig],'','w','HF Change','','','HorizontalAlignment','left','foregroundcolor',[1 1 1],'fontsize',9);
    k1=k1+1;l=l+ho2;
    h(k1)=add_button('text',[q1,l,q3,hig],'LF std:',[0.4 0.2 1],'LF STD text','','','foregroundcolor',[1 1 1],'fontsize',9,'visible','off');
    k1=k1+1;
    h(k1)=add_button('edit',[q4,l,q5,hig],'',[0.4 0.2 1],'LF STD Left','','','HorizontalAlignment','left','foregroundcolor',[1 1 1],'fontsize',9,'visible','off');
    k1=k1+1;
    h(k1)=add_button('edit',[q6,l,q5,hig],'',[0.4 0.2 1],'LF STD Right','','','HorizontalAlignment','left','foregroundcolor',[1 1 1],'fontsize',9,'visible','off');
    k1=k1+1;
    h(k1)=add_button('edit',[q7,l,q5,hig],'','w','LF STD Change','','','HorizontalAlignment','left','foregroundcolor',[1 1 1],'fontsize',9,'visible','off');
    k1=k1+1;l=l+ho2;
    h(k1)=add_button('text',[q1,l,q3,hig],'LF mean',[0.4 0.2 1],'','','','foregroundcolor',[1 1 1],'fontsize',9);
    k1=k1+1;
    h(k1)=add_button('edit',[q4,l,q5,hig],'',[0.4 0.2 1],'LF Left','','','HorizontalAlignment','left','foregroundcolor',[1 1 1],'fontsize',9);
    k1=k1+1;
    h(k1)=add_button('edit',[q6,l,q5,hig],'',[0.4 0.2 1],'LF Right','','','HorizontalAlignment','left','foregroundcolor',[1 1 1],'fontsize',9);
    k1=k1+1;
    h(k1)=add_button('edit',[q7,l,q5,hig],'','w','LF Change','','','HorizontalAlignment','left','foregroundcolor',[1 1 1],'fontsize',9);
    k1=k1+1;l=l+ho2;
    h(k1)=add_button('text',[q1,l,q3,hig],'Signal std:',[0.4 0.2 1],'Signal STD text','','','foregroundcolor',[1 1 1],'fontsize',9,'visible','off');
    k1=k1+1;
    h(k1)=add_button('edit',[q4,l,q5,hig],'',[0.4 0.2 1],'Signal STD Left','','','HorizontalAlignment','left','foregroundcolor',[1 1 1],'fontsize',9,'visible','off');
    k1=k1+1;
    h(k1)=add_button('edit',[q6,l,q5,hig],'',[0.4 0.2 1],'Signal STD Right','','','HorizontalAlignment','left','foregroundcolor',[1 1 1],'fontsize',9,'visible','off');
    k1=k1+1;
    h(k1)=add_button('edit',[q7,l,q5,hig],'','w','Signal STD Change','','','HorizontalAlignment','left','foregroundcolor',[1 1 1],'fontsize',9,'visible','off');
    k1=k1+1;l=l+ho2;
    h(k1)=add_button('text',[q1,l,q3,hig],'Signal mean:',[0.4 0.2 1],'','','','foregroundcolor',[1 1 1],'fontsize',9);
    k1=k1+1;
    h(k1)=add_button('edit',[q4,l,q5,hig],'',[0.4 0.2 1],'Signal Left','','','HorizontalAlignment','left','foregroundcolor',[1 1 1],'fontsize',9);
    k1=k1+1;
    h(k1)=add_button('edit',[q6,l,q5,hig],'',[0.4 0.2 1],'Signal Right','','','HorizontalAlignment','left','foregroundcolor',[1 1 1],'fontsize',9);
    k1=k1+1;
    h(k1)=add_button('edit',[q7,l,q5,hig],'','w','Signal Change','','','HorizontalAlignment','left','foregroundcolor',[1 1 1],'fontsize',9);
    k1=k1+1;l=l+ho2;

    h(k1)=add_button('text',[q1,l,q3,hig],'Diff [sec]',[0.4 0.2 1],'','','','foregroundcolor',[1 1 1],'fontsize',9);
    k1=k1+1;
    h(k1)=add_button('edit',[q4,l,q5,hig],'',[0.4 0.2 1],'Diff Left','','','HorizontalAlignment','left','foregroundcolor',[1 1 1],'fontsize',9);
    k1=k1+1;
    h(k1)=add_button('edit',[q6,l,q5,hig],'',[0.4 0.2 1],'Diff Right','','','HorizontalAlignment','left','foregroundcolor',[1 1 1],'fontsize',9);
    k1=k1+1;l=l+ho2;
    h(k1)=add_button('text',[q1,l,q3,hig],'End [sec]',[0.4 0.2 1],'','','','foregroundcolor',[1 1 1],'fontsize',9);
    k1=k1+1;
    h(k1)=add_button('edit',[q4,l,q5,hig],'',[0.4 0.2 1],'End Left','','','HorizontalAlignment','left','foregroundcolor',[1 1 1],'fontsize',9);
    k1=k1+1;
    h(k1)=add_button('edit',[q6,l,q5,hig],'',[0.4 0.2 1],'End Right','','','HorizontalAlignment','left','foregroundcolor',[1 1 1],'fontsize',9);
    k1=k1+1;l=l+ho2;
    h(k1)=add_button('text',[q1,l,q3,hig],'Start [sec]',[0.4 0.2 1],'','','','foregroundcolor',[1 1 1],'fontsize',9);
    k1=k1+1;
    h(k1)=add_button('edit',[q4,l,q5,hig],'',[0.4 0.2 1],'Start Left','','','HorizontalAlignment','left','foregroundcolor',[1 1 1],'fontsize',9);
    k1=k1+1;
    h(k1)=add_button('edit',[q6,l,q5,hig],'',[0.4 0.2 1],'Start Right','','','HorizontalAlignment','left','foregroundcolor',[1 1 1],'fontsize',9);
    k1=k1+1;l=l+ho2;

    
    h(k1)=add_button('text',[q1,l,q3,hig],'Region',[0.4 0.2 1],'','','','foregroundcolor',[1 1 1],'fontsize',9);
    k1=k1+1;
    h(k1)=add_button('text',[q4,l,q5,hig],'Left',[0.4 0.2 1],'','','','foregroundcolor',[1 1 1],'fontsize',9);
    k1=k1+1;
    h(k1)=add_button('text',[q6,l,q5,hig],'Right',[0.4 0.2 1],'','','','foregroundcolor',[1 1 1],'fontsize',9);
    k1=k1+1;
    h(k1)=add_button('text',[q7,l,q5,hig],'Change (%)',[0.4 0.2 1],'','','','foregroundcolor',[1 1 1],'fontsize',9);
    k1=k1+1;l=l+ho2;
    h(k1)=add_button('check',[q1,l,(q2-0.02)/3,hig],'Fix length:',[0.4 0.2 1],'Fix length check','','','foregroundcolor',[1 1 1],'fontsize',9,'enable','on','value',1);
    k1=k1+1;
    h(k1)=add_button('edit',[q1+(q2-0.02)/3+0.01,l,(q2-0.02)/3,hig],'300',[0.4 0.2 1],'Fix length edit','','','foregroundcolor',[1 1 1],'fontsize',9,'enable','on');
    k1=k1+1;
    h(k1)=add_button('radio',[q1+2*(q2-0.02)/3+0.02,l,(q2-0.02)/3,hig],'Equal regions',[0.4 0.2 1],'Equal Regions','','','foregroundcolor',[1 1 1],'fontsize',9,'value',1);
   k1=k1+1;l=l+ho2;
     h(k1)=add_button('push',[q1,l,(q2-0.02)/3,hig],'Select Regions',[0.4 0.2 1],'Select Regions','show_integrals_fcn(''Select Regions'');','','foregroundcolor',[1 1 1],'fontsize',9);
    k1=k1+1;
     h(k1)=add_button('push',[q1+(q2-0.02)/3+0.01,l,(q2-0.02)/3,hig],'Clear Regions',[0.4 0.2 1],'Clear Regions','show_integrals_fcn(''Clear Regions'');','','enable','off','foregroundcolor',[1 1 1],'fontsize',9);
    k1=k1+1;
    h(k1)=add_button('push',[q1+2*(q2-0.02)/3+0.02,l,(q2-0.02)/3,hig],'Show STD',[0.4 0.2 1],'Show STD','show_integrals_fcn(''Show STD'');','','enable','on','foregroundcolor',[1 1 1],'fontsize',9);
   a=add_frame(h);
    ap=get(a,'pos');

    k1=1;
    ho1=0.04;
    ho2=ho1+0.005;
    l=ap(2)+ap(4)+0.05;

    hh(k1)=add_button('check',[q1,l,q2,hig],'Log',[],'Log','show_integrals_fcn(''set Log'');');
    k1=k1+1;l=l+ho2;
    hh(k1)=add_button('radio',[q1,l,q2,hig],'X Zoom',[],'X Zoom','show_integrals_fcn(''set X Zoom'');');
    k1=k1+1;l=l+ho2;
    hh(k1)=add_button('radio',[q1,l,q2,hig],'Y Zoom',[],'Y Zoom','show_integrals_fcn(''set Y Zoom'');');
    k1=k1+1;l=l+ho2;
    hh(k1)=add_button('edit',[0.16,l,.05,hig],'35',[],'Median length','show_integrals_fcn(''update Median'');');
    k1=k1+1;
    hh(k1)=add_button('text',[q1,l,.14,hig],'Median Window:',[],'');
    k1=k1+1;l=l+ho2;
    hh(k1)=add_button('check',[q1,l,q2,hig],'Median Filter',[],'Median Filter','show_integrals_fcn(''set Median'');');
    set(hh,'backgroundcolor',[0.4 0.2 1],'foregroundcolor',[1 1 1],'fontsize',9); 
    show_integrals_fcn('set Y Zoom',h);
    set(findobj(gcf,'tag','Median length'),'background','w','ForegroundColor','k');
    add_frame(hh);

    add_icon([.01,.94,.04,.04],'print','printdlg(gcf);');
    add_icon([.055,.94,.04,.04],'help','helpdlg(help(''show_integrals''));');
    add_icon([.1,.94,.04,.04],'exit','show_integrals_fcn(''exit'');','tag','close');
    add_button('push',[.15,.94,.14,.04],'Export figure',[],'','export_fig(gcf);');
else
    varargout{1}=L;
    varargout{2}=H;
    varargout{3}=L./H;
end

function a=add_frame(h)
hh1=get(h,'pos');
hh2=cat(1,hh1{:});
m1=min(hh2(:,1));
m2=min(hh2(:,2));
m3=max(hh2(:,1)+hh2(:,3));
m4=max(hh2(:,2)+hh2(:,4));
a=add_button('frame',[m1-0.01,m2-0.01,m3-m1+0.02,m4-m2+0.02]);
set(a,'backgroundcolor',[.251,.502,.502],'foregroundcolor',[.251,.502,.502]);
change_order_of_object(a,'bottom');
    