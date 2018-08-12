function [x,y,ax] = get_points(varargin)

global GETPTS_FIG GETPTS_AX GETPTS_H1 GETPTS_H2 GETPTS_OBJ


if ((nargin >= 1) & (isstr(varargin{1})))
    feval(varargin{:});
    return;
end

if (nargin < 1)
    GETPTS_AX = [];
    GETPTS_FIG = gcf;
else
    if (~ishandle(varargin{1}))
        error('First argument is not a valid handle');
    end
    
    switch get(varargin{1}, 'Type')
        case 'figure'
            GETPTS_FIG = varargin{1};
            GETPTS_AX = [];
            GETPTS_OBJ='figure';
        case 'axes'
            GETPTS_AX = varargin{1};
            GETPTS_FIG = get(GETPTS_AX, 'Parent');
            GETPTS_OBJ='axes';      
        otherwise
            error('First argument should be a figure or axes handle');
    end
end

% Bring target figure forward
figure(GETPTS_FIG);

% Remember initial figure state
state = uisuspend(GETPTS_FIG);

% Set up initial callbacks for initial stage
pointerShape=[];
for k=1:length(varargin)
    if isstr(varargin{k})&strcmp(lower(varargin{k}),'pointer')
        pointerShape=varargin{k+1};
        varargin(k+(0:1))=[];
        set(GETPTS_FIG, 'WindowButtonDownFcn', 'get_points(''Button_Down'');', ...
            'KeyPressFcn', 'get_points(''KeyPress'');', ...
            'Pointer', pointerShape);            
        break
    end
end
if isempty(pointerShape)
    set(GETPTS_FIG, 'WindowButtonDownFcn', 'get_points(''Button_Down'');', ...
        'KeyPressFcn', 'get_points(''KeyPress'');');    
    set_pointer('getpoints');
end

set(GETPTS_FIG, 'WindowButtonMotionFcn','show_pointer_location(1);')

% Initialize the lines to be used for the drag
markerSize = 9;

GETPTS_H1=line('XData',[],'YData',[],'Visible','off','Clipping', 'off', ...
    'Color','c','LineStyle','none','Marker','+','MarkerSize',markerSize,'EraseMode', 'xor');

GETPTS_H2=line('XData',[],'YData',[],'Visible','off','Clipping', 'off', ...
    'Color','m','LineStyle','none','Marker','x','MarkerSize',markerSize,'EraseMode', 'xor');
if ~isempty(GETPTS_AX)
    set([GETPTS_H1,GETPTS_H2],'Parent',GETPTS_AX)
end
if length(varargin)==2 & ~ischar(varargin{2}(1))
    set(GETPTS_H1,'userdata',varargin{2})
end


% We're ready; wait for the user to do the drag
% Wrap the call to waitfor in try-catch so we'll
% have a chance to clean up after ourselves.
errCatch = 0;
try
    waitfor(GETPTS_H1, 'UserData', 'Completed');
catch
    errCatch=1;
end

% After the waitfor, if GETPTS_H1 is still valid
% and its UserData is 'Completed', then the user
% completed the drag.  If not, the user interrupted
% the action somehow, perhaps by a Ctrl-C in the
% command window or by closing the figure.

if (errCatch == 1)
    errStatus = 'trap';
    
elseif (~ishandle(GETPTS_H1) | ...
        ~strcmp(get(GETPTS_H1, 'UserData'), 'Completed'))
    errStatus = 'unknown';
    
else
    errStatus = 'ok';
    x = get(GETPTS_H1, 'XData');
    y = get(GETPTS_H1, 'YData');
    x = x(:);
    y = y(:);
    % If no points were selected, return rectangular empties.
    % This makes it easier to handle degenerate cases in
    % functions that call getpts.
    ax=[];
    if ~isempty(x) |~isempty(y)
        ax=get(GETPTS_H1,'parent');
    end
    if (isempty(x))
        x = zeros(0,1);
    end
    if (isempty(y))
        y = zeros(0,1);
    end
end

% Delete the animation objects
if (ishandle(GETPTS_H1))
    delete(GETPTS_H1);
end
if (ishandle(GETPTS_H2))
    delete(GETPTS_H2);
end

% Restore the figure state
if (ishandle(GETPTS_FIG))
    uirestore(state);
end

% Clean up the global workspace
clear global GETPTS_FIG GETPTS_AX GETPTS_H1 GETPTS_H2 GETPTS_OBJ
delete(findobj(gcf,'tag','pointLoc'))

% Depending on the error status, return the answer or generate
% an error message.
switch errStatus
    case 'ok'
        % No action needed.
        
    case 'trap'
        % An error was trapped during the waitfor
        error('Interruption during mouse point selection.');
        
    case 'unknown'
        % User did something to cause the point selection to
        % terminate abnormally.  For example, we would get here
        % if the user closed the figure in the middle of the selection.
        error('Interruption during mouse point selection.');
end


%--------------------------------------------------
% Subfunction KeyPress
%--------------------------------------------------
function KeyPress

global GETPTS_FIG GETPTS_AX GETPTS_H1 GETPTS_H2


key = real(get(GETPTS_FIG, 'CurrentCharacter'));
switch key
    case {8, 127}  % delete and backspace keys
        x = get(GETPTS_H1, 'XData');
        y = get(GETPTS_H1, 'YData');
        switch length(x)
            case 0
                % nothing to do
            case 1
                % remove point and start over
                set([GETPTS_H1 GETPTS_H2], ...
                    'XData', [], ...
                    'YData', []);
                set(GETPTS_FIG, 'WindowButtonDownFcn', ...
                    'get_points(''Button_Down'');');
            otherwise
                % remove last point
                set([GETPTS_H1 GETPTS_H2], ...
                    'XData', x(1:end-1), ...
                    'YData', y(1:end-1));
        end
        
    case {13, 3}   % enter and return keys
        % return control to line after waitfor
        set(GETPTS_H1, 'UserData', 'Completed');
        
end

%--------------------------------------------------
% Subfunction FirstButtonDown
%--------------------------------------------------
function Button_Down

global GETPTS_FIG GETPTS_AX GETPTS_H1 GETPTS_H2 GETPTS_OBJ

if isempty(GETPTS_AX)
    GETPTS_AX=get(GETPTS_FIG,'currentaxes');
end
set([GETPTS_H1,GETPTS_H2],'parent',GETPTS_AX)

[newx,newy] = get_current_point(GETPTS_AX);
if ~isempty(newx)
    if strcmp(GETPTS_OBJ,'axes')
        xl=get(GETPTS_AX,'xlim');
        yl=get(GETPTS_AX,'ylim');
        if newx<xl(1) |newx>xl(2) | newy<yl(1) |newy>yl(2) 
            newx=[];
        end
    end
end
x = get(GETPTS_H1, 'XData');
y = get(GETPTS_H2, 'YData');
if isempty(x)& isempty(newx),return,end
if isempty(x),x=newx;y=newy;else
    if ~isempty(newx),x=[x,newx];y=[y,newy];end
end


set([GETPTS_H1 GETPTS_H2], ...
    'XData', x, ...
    'YData', y, ...
    'Visible', 'on');

if (~strcmp(get(GETPTS_FIG, 'SelectionType'), 'normal'))
    % We're done!
    set(GETPTS_H1, 'UserData', 'Completed');
end
n=get(GETPTS_H1,'userdata');
if ~isempty(n)
    if length(x)==n
        set(GETPTS_H1, 'UserData', 'Completed');
    end
end


    %----------------------------------------------------
    % Subfunction Get_current_point
    %----------------------------------------------------
    function [x,y] = get_current_point(axHandle);
    if isempty(axHandle);
        a=findobj(gcbf,'type','axes');
        for k=1:length(a)
            Z = get(ax,'CurrentPoint');
            xlim = get(ax,'xlim');
            ylim = get(ax,'ylim');
            if (xlim(1) <= Z(1,1) & Z(1,1) <= xlim(2) &  ylim(1) <= Z(1,2) & Z(1,2) <= ylim(2))
                break
            end
        end
        axHandle=a(k); 
    end
    pt = get(axHandle, 'CurrentPoint');
    x = pt(1,1);
    y = pt(1,2);
    
    % What is the extent of the idealized screen pixel in axes
    % data space?
    
    axUnits = get(axHandle, 'Units');
    set(axHandle, 'Units', 'pixels');
    axPos = get(axHandle, 'Position');
    set(axHandle, 'Units', axUnits);
    
    axPixelWidth = axPos(3);
    axPixelHeight = axPos(4);
    
    axXLim = get(axHandle, 'XLim');
    axYLim = get(axHandle, 'YLim');
    
    xExtentPerPixel = abs(diff(axXLim)) / axPixelWidth;
    yExtentPerPixel = abs(diff(axYLim)) / axPixelHeight;
    
    x = x + xExtentPerPixel/2;
    y = y + yExtentPerPixel/2;