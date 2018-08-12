function a=add_icon(pos,icon,cb,varargin)
% This function draws a pushbutton with image on it (active icon).
% You must provide the position and the type for the icon. Other input parameters are optional.
% Icon data is a mat file, located in SPAT\ICONS (you need to provide the icon file name. 

a=add_button('pushbutton',pos);
if nargin>2
    set(a,'callback',cb)
end
icons_folder=get_icons_folder;

load([icons_folder,lower(icon),'_icon.mat'])
set(a,'cdata',x)
for k=1:2:length(varargin)
    set(a,varargin{k},varargin{k+1})
end

function icons_folder=get_icons_folder;
c=which('spat.m');
if isempty(c),
    c=which('help_icon.mat');
else
    icons_folder=[get_file_prop(c,'path'),'icons\'];
    return
end
if isempty(c), c=which('help_icon.mat'),;end
if isempty(c), c=which('psd_fcn.m'),;end
if isempty(c), c=which('add_button.m'),;end
D=cd;
a=get_file_prop(c,'path');
cd(a)
cd('..');
c=cd;
cd(D)
icons_folder=[c,'\icons\'];