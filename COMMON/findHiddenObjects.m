function obj=findHiddenObjects(varargin)
if ishandle(varargin{1}),
    obj=varargin{1};
    varargin(1)=[];
else
    obj=0;
end
set(0,'showhiddenhandles','on')
obj=findobj(obj,varargin{1:2:end},varargin{2:2:end});
for k=1:2:length(varargin)
    obj=findobj(obj,varargin{k},varargin{k+1});
end
set(0,'showhiddenhandles','off')
