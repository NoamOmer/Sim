function SelectColor(tag)
S=get(gcbo,'string');
v=get(gcbo,'value');
Col=S{v};
switch lower(Col)
    case {'red','r'}
        col=[1,0,0];
    case {'green','g'}
        col=[0,1,0];
    case {'blue','b'}
        col=[0,0,1];
    case {'white','w'}
        col=[1,1,1];
    case {'black','k'}
        col=[0,0,0];
    case {'yellow','y'}
        col=[0,1,1];
    case {'none','n'}
        col=[];
    case {'custom','customize'}
        col=uisetcolor;
end
obj=findobj(gcbf,'tag',tag);
set(obj,'userdata',col)
if isempty(col)
    X='none';
else
    X='';
    for k=1:length(col)
        X=[X,',',num2str(col(k),3)];
    end
    X(1)=[];
end
if strcmp(lower(get(obj,'style')),'edit')
    set(obj,'string',X)
end
if strcmp(lower(get(obj,'style')),'text')
    set(obj,'string',X)
end