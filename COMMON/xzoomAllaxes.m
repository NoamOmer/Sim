function out = xzoomallaxes(varargin)
fig=varargin{1};
if ~ishandle(fig)
    error('First argument must be a figure handle.')
end
zoomCommand=lower(varargin{2});

%
% zoomCommand 'off'
%
if strcmp(zoomCommand,'off')
    % turn zoom off
    doZoomOff(fig);
    scribefiglisten(fig,'off');
    state = getappdata(fig,'ZOOMFigureState');
    if ~isempty(state)
        ptr = get(fig,'pointer');
        uirestore(state,'nouicontrols');
        set(fig,'pointer',ptr)
        if isappdata(fig,'ZOOMFigureState')
            rmappdata(fig,'ZOOMFigureState');
        end
        if isappdata(fig,'ZoomOnState')
            rmappdata(fig,'ZoomOnState');
        end
        
    end
    % done, go home.
    return
end

% set some things we need for other zoomCommands
ax=get(fig,'currentaxes');
rbbox_mode = 0;
% initialize unconstrained state
zoomx = 1; zoomy = 1;
%
% the zoomCommand is 'toggle'
%
if strcmp(zoomCommand,'toggle'),
    state = getappdata(fig,'ZOOMFigureState');
    if isempty(state)
        xzoomallaxes(fig,'on');
    else
        xzoomallaxes(fig,'off');
    end
    % done, go home
    return
end % if


% Set zoomx,zoomy and zoomCommand for constrained zooms
if strcmp(zoomCommand,'xdown'),
    zoomy = 0; zoomCommand = 'down'; % Constrain y
elseif strcmp(zoomCommand,'ydown')
    zoomx = 0; zoomCommand = 'down'; % Constrain x
end

switch zoomCommand
case 'down'
    % Activate axis that is clicked in
    allAxes = findall(datachildren(fig),'flat','type','axes');
    ZOOM_found = 0;
    
    % this test may be causing failures for 3d axes
    for i=1:length(allAxes)
        ax=allAxes(i);
        ZOOM_Pt1 = get(ax,'CurrentPoint');
        xlim = get(ax,'xlim');
        ylim = get(ax,'ylim');
        if (xlim(1) <= ZOOM_Pt1(1,1) & ZOOM_Pt1(1,1) <= xlim(2) & ...
                ylim(1) <= ZOOM_Pt1(1,2) & ZOOM_Pt1(1,2) <= ylim(2))
            ZOOM_found = 1;
            set(fig,'currentaxes',ax);
            break
        end
    end
    
    if ZOOM_found==0
        return
    end
    
    % Check for selection type
    selection_type = get(fig,'SelectionType');
    zoomMode = getappdata(fig,'ZOOMFigureMode');
    
    axz = get(ax,'ZLabel');
    
    
    if isempty(zoomMode) | strcmp(zoomMode,'in');
        switch selection_type
        case 'normal'
            % Zoom in
            m = 1;
            scale_factor = 2; % the default zooming factor
        case 'open'
            % Zoom all the way out
            xzoomallaxes(fig,'out');
            return;
        otherwise
            % Zoom partially out
            m = -1;
            scale_factor = 2;
        end
    elseif strcmp(zoomMode,'out')
        switch selection_type
        case 'normal'
            % Zoom partially out
            m = -1;
            scale_factor = 2;
        case 'open'
            % Zoom all the way out
            xzoomallaxes(fig,'out');
            return;
        otherwise
            % Zoom in
            m = 1;
            scale_factor = 2; % the default zooming factor
        end
    else % unrecognized zoomMode
        return
    end
    
    ZOOM_Pt1 = get_currentpoint(ax);
    ZOOM_Pt2 = ZOOM_Pt1;
    center = ZOOM_Pt1;
    
    if (m == 1)
        % Zoom in
        units = get(fig,'units'); set(fig,'units','pixels')
        
        rbbox([get(fig,'currentpoint') 0 0],get(fig,'currentpoint'),fig);
        
        ZOOM_Pt2 = get_currentpoint(ax);
        set(fig,'units',units)
        
        % Note the currentpoint is set by having a non-trivial up function.
        if min(abs(ZOOM_Pt1-ZOOM_Pt2)) >= ...
                min(.01*[diff(get_xlim(ax)) diff(get_ylim(ax))]),
            % determine axis from rbbox 
            a = [ZOOM_Pt1;ZOOM_Pt2]; a = [min(a);max(a)];
            
            % Undo the effect of get_currentpoint for log axes
            if strcmp(get(ax,'XScale'),'log'),
                a(1:2) = 10.^a(1:2);
            end
            if strcmp(get(ax,'YScale'),'log'),
                a(3:4) = 10.^a(3:4);
            end
            rbbox_mode = 1;
        end
    end
    limits = xzoomallaxes(fig,'getlimits');
case 'on',
    state = getappdata(fig,'ZOOMFigureState');
    if isempty(state),
        % turn off all other interactive modes
        state = uiclearmode(fig,'docontext','zoom',fig,'off');
        % restore button down functions for uicontrol children of the figure
        uirestore(state,'uicontrols');
        setappdata(fig,'ZOOMFigureState',state);
    end
    
    set(fig,'windowbuttondownfcn','xzoomAllaxes(gcbf,''down'')', ...
        'windowbuttonupfcn','ones;', ...
        'windowbuttonmotionfcn','', ...
        'buttondownfcn','', ...
        'interruptible','on');
    set(ax,'interruptible','on');
    % set an appdata so it will always be possible to 
    % determine whether or not zoom is on and in what
    % type of 'on' state it is.
    % this appdata will not exist when zoom is off
    setappdata(fig,'ZoomOnState','on');
    
    scribefiglisten(fig,'on');
    
    doZoomIn(fig)
    return
        
case 'getmode'
    state = getappdata(fig,'ZOOMFigureState');
    if isempty(state)
        out = 'off';
    else
        mode = getappdata(fig,'ZOOMFigureMode');
        if isempty(mode)
            out = 'xon';
        else
            out = mode;
        end
    end
    return

case 'inmode'
    xzoomAllaxes(fig,'on');
    doZoomIn(fig)
    return   
  
case 'reset',
    axz = get(ax,'ZLabel');
    if isappdata(axz,'ZOOMAxesData')
        rmappdata(axz,'ZOOMAxesData');
    end
    return
    
case 'out',
    limits = xzoomAllaxes(fig,'getlimits');
    center = [sum(get_xlim(ax))/2 sum(get_ylim(ax))/2];
    m = -inf; % Zoom totally out
    
case 'getlimits', % Get axis limits
      axz = get(ax,'ZLabel');
        limits = getappdata(axz,'ZOOMAxesData');
        % Do simple checking of userdata
        if size(limits,2)==4 & size(limits,1)<=2, 
            if all(limits(1,[1 3])<limits(1,[2 4])), 
                getlimits = 0; out = limits(1,:);
                return   % Quick return
            else
                getlimits = -1; % Don't munge data
            end
        else
            if isempty(limits)
                getlimits = 1;
            else 
                getlimits = -1;
            end
        end
  
        % If I've made it to here, we need to compute appropriate axis
        % limits.
  
        if isempty(getappdata(axz,'ZOOMAxesData')),
            % Use quick method if possible
            xlim = get_xlim(ax); xmin = xlim(1); xmax = xlim(2); 
            ylim = get_ylim(ax); ymin = ylim(1); ymax = ylim(2); 
      
        elseif strcmp(get(ax,'xLimMode'),'auto') & ...
                strcmp(get(ax,'yLimMode'),'auto'),
            % Use automatic limits if possible
            xlim = get_xlim(ax); xmin = xlim(1); xmax = xlim(2); 
            ylim = get_ylim(ax); ymin = ylim(1); ymax = ylim(2); 
      
        else
            % Use slow method only if someone else is using the userdata
            h = get(ax,'Children');
            xmin = inf; xmax = -inf; ymin = inf; ymax = -inf;
            for i=1:length(h),
                t = get(h(i),'Type');
                if ~strcmp(t,'text'),
                    if strcmp(t,'image'), % Determine axis limits for image
                        x = get(h(i),'Xdata'); y = get(h(i),'Ydata');
                        x = [min(min(x)) max(max(x))];
                        y = [min(min(y)) max(max(y))];
                        [ma,na] = size(get(h(i),'Cdata'));
                        if na>1 
                            dx = diff(x)/(na-1);
                        else 
                            dx = 1;
                        end
                        if ma>1
                            dy = diff(y)/(ma-1);
                        else
                            dy = 1;
                        end
                        x = x + [-dx dx]/2; y = y + [-dy dy]/2;
                    end
                    xmin = min(xmin,min(min(x)));
                    xmax = max(xmax,max(max(x)));
                    ymin = min(ymin,min(min(y)));
                    ymax = max(ymax,max(max(y)));
                end
            end
      
            % Use automatic limits if in use (override previous calculation)
            if strcmp(get(ax,'xLimMode'),'auto'),
                xlim = get_xlim(ax); xmin = xlim(1); xmax = xlim(2); 
            end
            if strcmp(get(ax,'yLimMode'),'auto'),
                ylim = get_ylim(ax); ymin = ylim(1); ymax = ylim(2); 
            end
        end
  
        limits = [xmin xmax ymin ymax];
        if getlimits~=-1, % Don't munge existing data.
                            % Store limits ZOOMAxesData
                            % store it with the ZLabel, so that it's cleared if the 
                            % user plots again into this axis.  If that happens, this
                            % state is cleared
                            axz = get(ax,'ZLabel');
                            setappdata(axz,'ZOOMAxesData',limits);
        end
  
        out = limits;
        return
  otherwise
    error(['Unknown option: ',zoomCommand,'.']);
end

%---------------------------------------------------------------------------------%

%
% Actual zoom operation
%

if ~rbbox_mode,
    xmin = limits(1); xmax = limits(2); 
    ymin = limits(3); ymax = limits(4);
    
    if m==(-inf),
        dx = xmax-xmin;
        dy = ymax-ymin;
    else
        dx = diff(get_xlim(ax))*(scale_factor.^(-m-1)); dx = min(dx,xmax-xmin);
        dy = diff(get_ylim(ax))*(scale_factor.^(-m-1)); dy = min(dy,ymax-ymin);
    end
    
    % Limit zoom.
    center = max(center,[xmin ymin] + [dx dy]);
    center = min(center,[xmax ymax] - [dx dy]);
    a = [max(xmin,center(1)-dx) min(xmax,center(1)+dx) ...
            max(ymin,center(2)-dy) min(ymax,center(2)+dy)];
    
    % Check for log axes and return to linear values.
    if strcmp(get(ax,'XScale'),'log'),
        a(1:2) = 10.^a(1:2);
    end
    if strcmp(get(ax,'YScale'),'log'),
        a(3:4) = 10.^a(3:4);
    end
    
end

% Check for axis equal and update a as necessary
if strcmp(get(ax,'plotboxaspectratiomode'),'manual') & ...
        strcmp(get(ax,'dataaspectratiomode'),'manual')
    ratio = get(ax,'plotboxaspectratio')./get(ax,'dataaspectratio');
    dx = a(2)-a(1);
    dy = a(4)-a(3);
    [kmax,k] = max([dx dy]./ratio(1:2));
    if k==1
        dy = kmax*ratio(2);
        a(3:4) = mean(a(3:4))+[-dy dy]/2;
    else
        dx = kmax*ratio(1);
        a(1:2) = mean(a(1:2))+[-dx dx]/2;
    end
end
if a(1)==a(2), return, end 
allAxes=findobj(fig,'type','axes');
set(allAxes,'ylimmode','manual','xlimmode','manual')
subMenu=findHiddenObjects(fig,'tag','subMenu');
opt=get(findHiddenObjects(subMenu,'state','on'),'tag');
switch opt
case 'zoomGcaXY'
    set(ax,'xlim',a(1:2),'ylim',a(3:4));
case 'zoomGcaX'
    set(ax,'xlim',a(1:2));
case 'zoomGcaY'
    set(ax,'ylim',a(3:4));
case 'zoomXYX'
    set(ax,'ylim',a(3:4));
    set(allAxes,'xlim',a(1:2));
case 'zoomXX'
    set(allAxes,'xlim',a(1:2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = get_currentpoint(ax)
%GET_CURRENTPOINT Return equivalent linear scale current point
p = get(ax,'currentpoint'); p = p(1,1:2);
if strcmp(get(ax,'XScale'),'log'),
    p(1) = log10(p(1));
end
if strcmp(get(ax,'YScale'),'log'),
    p(2) = log10(p(2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xlim = get_xlim(ax)
%GET_XLIM Return equivalent linear scale xlim
xlim = get(ax,'xlim');
if strcmp(get(ax,'XScale'),'log'),
    xlim = log10(xlim);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ylim = get_ylim(ax)
%GET_YLIM Return equivalent linear scale ylim
ylim = get(ax,'ylim');
if strcmp(get(ax,'YScale'),'log'),
    ylim = log10(ylim);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doZoomIn(fig)
set(findall(fig,'Tag','xZoomInAll'),'State','on');   
setappdata(fig,'ZOOMFigureMode','in');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doZoomOff(fig)
if ~isempty(getappdata(fig,'ZoomFigureMode'))
    rmappdata(fig,'ZOOMFigureMode');
end