function a=add_button(style,pos,str,color,tag,call,userdata,varargin);
% This function draws a uicontrol on the current axis.
% You must provide the style and position for the button. Other input parameters are optional.

% Define default color
if findstr(lower(style),'edit') 
    but_color=[1 1 1];%White
elseif findstr(lower(style),'popup')
    but_color=[1 1 1];%white
else
    but_color=[.7,.7,.7];%Gray
end

if nargin<3
   str='';color=but_color;tag='';call='';userdata=[];
elseif nargin<4
   color=but_color;tag='';call='';userdata=[];
elseif nargin<5
   tag='';call='';userdata=[];  
elseif nargin<6
   call='';userdata=[];
elseif nargin<7
   userdata=[];
end
if isempty(color),color=but_color;end
if ismatrix(color)
    bgcolor=color(1,:);
    fgcolor=color(2,:);
else
    bgcolor=color;
    fgcolor=[0 0 0];    
end
a=uicontrol('style',style,'units','normalized','position',pos,'string',str,...
   'backgroundcolor',bgcolor,'tag',tag,'callback',call,'userdata',userdata,'foregroundcolor',fgcolor);

if strcmp(lower(style),'frame')
   set(a,'foregroundcolor',color);
end
for i=1:2:length(varargin)
   set(a,varargin{i},varargin{i+1})
end
 set(a,'tooltip',call(:)')
