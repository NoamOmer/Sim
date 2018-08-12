function PatchEditor(obj)
if ishandle(obj)
    f=createOptFigure(obj);
    return
end
opt=obj;
fig=gcbf;
switch opt
    case 'facecolor'
        SelectColor('FaceColor');
    case 'edgecolor'
        SelectColor('EdgeColor');
    case 'continue'
        obj=get(fig,'userdata');
        [pattern,mark]=getSelectedPattern(fig);
        faceColor=get(findobj(fig,'tag','FaceColor'),'userdata');   
        EdgeColor=get(findobj(fig,'tag','EdgeColor'),'userdata');
        lineWidth=str2num(get(findobj(fig,'tag','EdgeWidth'),'string'));
        delete(fig);
        a=get(obj,'userdata');
        if ishandle(a),delete(a),end
        if ~isempty(pattern)
            a=changePatchpattern(obj,pattern,faceColor,mark);
            faceColor='none';
        end
        if isempty(EdgeColor)
            EdgeColor='none';
        end
        if isempty(faceColor)
            faceColor='none';
        end
        set(obj,'userdata',a,'facecolor',faceColor,'edgecolor',EdgeColor,'linewidth',lineWidth)
        set(0,'showhiddenhandles','off')
        return        
    case 'cancel'
        delete(fig)        
        set(0,'showhiddenhandles','off')
        return
end
%**************************************************************************
%**************************************************************************
function fig=createOptFigure(obj);
N=32;
X=vec2mat(1:N*N,N,0);
X1=fliplr(X);
K=1:N;K=K(:);
A=[];B=[];C=[];D=[];E=[];F=[];
for k=-(N-1):6:(N-1)
    A{end+1}=diag(X,k);
    B{end+1}=diag(X1,k);
end
E={X(3:6:N,:)';X1(:,3:6:N)};
F={X(3:10:N,:)';X1(:,3:10:N)};
G={X(1:5:N,1:5:N)'};
H={X(2:7:N,2:7:N)'};

xl=0.5:0.1:1;
d=0.09;
fig=figure('windowstyle','normal','menubar','none','units','normalized','position',...
    [0.35,0.3,0.4,0.5],'userdata',obj);
add_button('text',[0.02,0.8,0.37,d],'Select surface patten');
a(1)=add_button('toggle',[0.4,0.8,d,d],'',[],'pattern','rad_tog',{X});
a(2)=add_button('toggle',[xl(1),0.7,d,d],'FullNegLines',[],'pattern','rad_tog',A);
a(3)=add_button('toggle',[xl(1),0.8,d,d],'NegLines',[],'pattern','rad_tog',A(1:2:end));
a(4)=add_button('toggle',[xl(2),0.7,d,d],'FullPosLines',[],'pattern','rad_tog',B);
a(5)=add_button('toggle',[xl(2),0.8,d,d],'PosLines',[],'pattern','rad_tog',B(1:2:end));
a(6)=add_button('toggle',[xl(3),0.7,d,d],'FullCrossLines',[],'pattern','rad_tog',[A,B]);
a(7)=add_button('toggle',[xl(3),0.8,d,d],'CrossLines',[],'pattern','rad_tog',[A(1:2:end),B(1:2:end)]);
a(8)=add_button('toggle',[xl(4),0.7,d,d],'FullSquareLines',[],'pattern','rad_tog',E);
a(9)=add_button('toggle',[xl(4),0.8,d,d],'SquareLines',[],'pattern','rad_tog',F);
a(10)=add_button('toggle',[xl(5),0.7,d,d],'FullDots',[],'pattern','rad_tog',G);
a(11)=add_button('toggle',[xl(5),0.8,d,d],'Dots',[],'pattern','rad_tog',H);
set(a,'ForegroundColor',get(a(1),'backgroundcolor'))
set(a(1),'value',1,'enable','inactive');
d=0.06;
color={'None','Red','Blue','Green','White','Black','Yellow','Customize'};
% Face color
[faceColor,fc,k]=get_current_color(obj,'facecolor');
add_button('text',[0.02,0.55,0.37,d],'Select surface color');
add_button('popup',[0.4,0.55,0.3,d],color,[],'','patcheditor(''facecolor'');',[],'value',k);
add_button('edit',[0.75,0.55,0.2,d],faceColor,[1 1 1],'FaceColor','',fc);
% Edge color
[edgeColor,ec,k]=get_current_color(obj,'edgecolor');
add_button('text',[0.02,0.45,0.37,d],'Select edge color');
add_button('popup',[0.4,0.45,0.3,d],color,[],'','patcheditor(''edgecolor'');',[],'value',k);
add_button('edit',[0.75,0.45,0.2,d],edgeColor,[1 1 1],'EdgeColor','',ec);
% Edge width
w=get(obj,'linewidth');
add_button('text',[0.02,0.35,0.37,d],'Select edge width');
add_button('edit',[0.4,0.35,0.1,d],num2str(w),[1 1 1],'EdgeWidth');
add_button('text',[0.51,0.35,0.2,d],'points');

add_button('push',[0.2,0.01,0.28,0.06],'Continue',[],'','patcheditor(''continue'');');
add_button('push',[0.52,0.01,0.28,0.06],'Cancel',[],'','patcheditor(''cancel'');');
refreshParameters(fig);
%**************************************************************************
%**************************************************************************
function refreshParameters(fig)
N=32;
X=vec2mat(1:N*N,N,0);
edge=[X(:,1);X(1,:)';X(:,end);X(end,:)'];
a=findobj(fig,'tag','pattern');

FC=[1 0 0];%get(findobj(fig,'tag','FaceColor'),'userdata');
EC=[0 0 1];%get(findobj(fig,'tag','EdgeColor'),'userdata');
for k=1:length(a)
    R=NaN*ones(N);G=R;B=R;
    x=get(a(k),'userdata');
    X=[];
    for m=1:length(x)
        s=x{m};
        X=[X;s(:)];
    end
    R(X)=FC(1);R(edge)=EC(1);
    G(X)=FC(2);G(edge)=EC(2);
    B(X)=FC(3);B(edge)=EC(3);
    C(:,:,1)=R;    C(:,:,2)=G;    C(:,:,3)=B;
    set(a(k),'cdata',C);
end
%**************************************************************************
%**************************************************************************
function [COL,col,k]=get_current_color(obj,prop)
col=get(obj,prop);
if isstr(col),
    COL='';k=1;return
else
    COL=num2str(col);
end
c=[NaN,NaN,NaN
    1,0,0
    0,0,1
    0,1,0
    1,1,1
    0,0,0
    1,1,0];
k=find(c(:,1)==col(1) & c(:,2)==col(2) & c(:,3)==col(3));
if isempty(k),k=8;end

function [pattern,mark]=getSelectedPattern(fig);
a=findobj(fig,'tag','pattern','value',1);
pattern=get(a,'string');
mark='.';