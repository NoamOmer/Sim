function putdowntext(action,fig)

if isempty(fig) | ~ishandle(fig), return, end   
toolButton = getappdata(fig,'ScribeCurrentToolButton');

if ~isempty(toolButton)  % aborting one operation and starting another
    if ishandle(toolButton) ...
            & strcmp(get(toolButton,'Type'),'uitoggletool') ...
            & toolButton ~= gcbo  % not the same button
        set(toolButton,'State','off');
    else
        toolButton = [];
    end
end
toolButton = action;

setappdata(fig,'ScribeCurrentToolButton',toolButton);

switch action
case 'zoomin'  
    switch xzoomAllaxes(fig,'getmode')
    case {'off' 'out'}
        xzoomAllaxes(fig,'inmode');
    case {'in' 'on'}
        xzoomAllaxes(fig,'off');
    end    
case 'zoomout'    
    switch xzoomAllaxes(fig,'getmode')
    case {'on' 'off' 'in'}
        xzoomAllaxes(fig,'outmode');
    case 'out'
        xzoomAllaxes(fig,'off');
    end
end