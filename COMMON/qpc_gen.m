function [Z,PHASE]=qpc_gen(varargin)
% Optional input parameters:
% Nqpc: # of qpc triads
% QPC matrices of [amplitudes(:),frequencies(:)]
% Nqpc+1) : uncoupled triads
% nlen : length of realizations
% Nr: # of realizations
% SF: desired sampling frequency
% Anoise: amplitude of noise
% Z=qpc_gen(2,[1,.1;1,.2;1,.3],[1,.55;1,.87;1,.1.42],[],512,50,10,0);

Z=[];
switch nargin
    case 0  
        [Nqpc,Ac,Fc,Nuc,Auc,Fuc,nlen,Nr,SF,Anoise]=add_figure;
    case 1
        feval(varargin{1});
        return
    otherwise   
        Nqpc=varargin{1};
        Ac=[];Fc=[];
        if Nqpc>0
            for k=1:Nqpc
                X=varargin{k+1};
                Ac(:,k)=X(:,1);
                Fc(:,k)=X(:,2);
            end
            Fc(3,:)=Fc(1,:)+Fc(2,:);
        end
        m=2;
        X=varargin{Nqpc+m};
        if isempty(X);Auc=[];Fuc=[];
        else,   Auc=X(:,1);   Fuc=X(:,2);end
        Nuc=size(Auc,1);
        l=Nqpc+m+1;
        nlen=varargin{l};
        Nr=varargin{l+1};
        SF=varargin{l+2};
        Anoise=varargin{l+3};
end

tvec = (0:nlen-1)'/SF; 
alpha = ones(nlen,1);
freqs=2*pi*[Fc(:);Fuc(:)]';
amps=[Ac(:);Auc(:)];
%[[Fc(:);Fuc(:)],[Ac(:);Auc(:)]]


Z=zeros(nlen,Nr);  
thq=[];
thu=[];
for k=1:Nr
    if Nqpc>0
        thq=(rand(Nqpc*3,1)-.5)*2*pi;
        thq=reshape(thq,Nqpc,3);
        %         thq(:,3) = thq(:,1) + thq(:,2);
        zz(k)=(rand(1,Nqpc)-.5)*pi*0;
        thq(:,3) = thq(:,1) + thq(:,2)+zz(k);
        thq=thq';
    end
    if Nuc>0
        thu=(rand(Nuc,1)-.5)*2*pi;
    end
    theta = [thq(:); thu(:)]';
    PHASE(:,k)=theta';
    y=zeros(nlen,1);
    if ~isempty(theta)
        y = cos(tvec*freqs+alpha*theta)*amps;
    end
    if Anoise>0
        y=y+Anoise*randn(nlen,1);
    end   
    Z(:,k) = y;
end
%*********************************************************************************
%
%*********************************************************************************
function [Nqpc,Ac,Fc,Nuc,Auc,Fuc,nlen,Nr,SF,Anoise]=add_figure;
fig=figure('name','QPC','windowstyle','modal');
add_button('frame',[.01,.77,.32,.19],'',[.7,0,1]);
add_button('text',[.02,.9,.2,.05],'N realizations');
N_R=add_button('edit',[.23,.9,.09,.05],'100',[1 1 1],'Nr');
add_button('text',[.02,.84,.2,.05],'Samples/realization');
NLEN=add_button('edit',[.23,.84,.09,.05],'1024',[1 1 1],'len');
add_button('text',[.02,.78,.2,.05],'Sampling Frequency');
S_F=add_button('edit',[.23,.78,.09,.05],'10',[1 1 1],'SF','qpc_gen(''change_sf'');',10);
call='qpc_gen(''ena_dis'');';
A=ones(4,3);
SF=10;
F(:,1)=[.02;.05;.15;.22];
F(:,2)=[.041;.044;.12;.12];
F(:,3)=F(:,1)+F(:,2);
F=F*SF;
%add_button(
for k=1:4
    STR=['qpc',num2str(k)];ASTR=['Aqpc',num2str(k)];
    dx=.09;K=.8-k*.14;
    qpc=['S1=findobj(gcf,''tag'',''',STR,''',''style'',''edit'');'...
            'S2=findobj(gcf,''tag'',''',STR,''',''style'',''text'');'...
            'F3=num2str(str2num(get(S1(1),''string''))+str2num(get(S1(2),''string'')));'...
            'set(S2,''string'',F3)'];
    add_button('frame',[.01,K-0.07,.51,.13],'',[.7,0,1]);
    cb_qpc(k)=add_button('checkbox',[.02,K,dx,.05],['QPC',num2str(k)],[0 0 1],'',call,...
        STR,'value',1);
    add_button('text',[.12,K,dx,.05],'Freq:',[.7 .7 .7]);
    freq(k,1)=add_button('edit',[.22,K,dx,.05],F(k,1),[1 1 1],STR,qpc);
    freq(k,2)=add_button('edit',[.32,K,dx,.05],F(k,2),[1 1 1],STR,qpc);
    freq(k,3)=add_button('text',[.42,K,dx,.05],F(k,3),[1 1 1],STR);
    K=K-.06;
    add_button('text',[.12,K,dx,.05],'Amp :',[.7 .7 .7]);
    ampqpc(k,1)=add_button('edit',[.22,K,dx,.05],A(k,1),[1 1 1],ASTR);
    ampqpc(k,2)=add_button('edit',[.32,K,dx,.05],A(k,2),[1 1 1],ASTR);
    ampqpc(k,3)=add_button('edit',[.42,K,dx,.05],A(k,3),[1 1 1],ASTR);
end
F=[.01;.03;.18;.3]*SF;
for k=1:4
    STR=['uc',num2str(k)];ASTR=['Auc',num2str(k)];
    dx=.09;K=.8-k*.14;
    add_button('frame',[.6,K-0.07,.31,.13],'',[.7,0,1]);
    cb_uc(k)=add_button('checkbox',[.61,K,dx,.05],['UC',num2str(k)],[],'',call,STR,...
        'value',0);
    add_button('text',[.71,K,dx,.05],'Freq:',[.7 .7 .7]);
    freq_uc(k)=add_button('edit',[.81,K,dx,.05],F(k),[.8 .8 .8],STR,qpc,'','enable','off');
    K=K-.06;
    add_button('text',[.71,K,dx,.05],'Amp :',[.7 .7 .7]);
    amp_uc(k)=add_button('edit',[.81,K,dx,.05],A(k,1),[.8 .8 .8],ASTR,'',[],'enable','off');
end

add_button('frame',[.01,K-.11,.33,.07],'',[.7,0,1]);
add_button('text',[.02,K-.1,.2,.05],'White noise amp:');
anoise=add_button('edit',[.23,K-.1,.1,.05],0,[1 1 1],'noise');
completed=add_button('pushbutton',[.82,.01,.16,.05],'Generate',[0 1 1],'',...
    'set(gco,''userdata'',''completed'');');
k=findobj(fig,'style','checkbox');
set(k,'backgroundcolor',[0 0 1],'foregroundcolor',[1 1 1])
k=findobj(fig,'style','text','tag','');
set(k,'backgroundcolor',[0 .6 1],'foregroundcolor',[1 1 1])


waitfor(completed,'userdata','completed');


if isempty(fig),return;end
Fc=[];Fuc=[];
Ac=[];Auc=[];
nlen=str2num(get(NLEN,'string'));
Nr=str2num(get(N_R,'string'));
SF=str2num(get(S_F,'string'));
l=0;
for k=1:4
    if get(cb_qpc(k),'value')==1
        l=l+1;F=[];
        for m=1:3
            F(m)=str2num(get(freq(k,m),'string'));
            AA(m)=str2num(get(ampqpc(k,m),'string'));
        end 
        [ll,I]=sort(F);
        Fc(l,:)=ll;
        Ac(l,:)=AA(I);
    end   
end
Fc=Fc';Ac=Ac';
Nqpc=l;
%Fc(:,3)=Fc(:,1)+Fc(:,2);
Nuc=4;
l=0;
for k=1:Nuc
    if get(cb_uc(k),'value')==1      
        l=l+1;
        Fuc(l)=str2num(get(freq_uc(k),'string'));
        Auc(l)=str2num(get(amp_uc(k),'string'));   
    end
end
Anoise=str2num(get(anoise,'string'));
Nuc=l;
close(fig);

%*********************************************************************************
%
%*********************************************************************************
function change_sf
oSF=get(gcbo,'userdata');
nSF=str2num(get(gcbo,'string'));
set(gcbo,'userdata',nSF);
for k=1:4
    S=findobj(gcf,'tag',['qpc',num2str(k)]);
    for l=1:3
        f=str2num(get(S(l),'string'))*nSF/oSF;
        set(S(l),'string',f);
    end
end

%*********************************************************************************
%
%*********************************************************************************
function ena_dis
T=lower(get(gcbo,'string'));
h=findobj(gcf,'tag',T);h1=findobj(gcf,'tag',['A',T]);
h=[h(:);h1(:)];
if get(gcbo,'value')==1,
    set(h,'enable','on','backgroundcolor',[1 1 1]) 
else,
    set(h,'enable','off','backgroundcolor',[.8,.8,.8])
end



