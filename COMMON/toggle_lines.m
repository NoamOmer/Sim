function toggle_lines(x,opt,fig)
if nargin<3
    fig=gcbf;
end
if isstr(x)
    x=findobj(fig,'tag',x);
end
switch lower(opt)
    case 'visible'        
        s=get(x(1),'visible');
        switch lower(s)
            case 'on'
                V='off';
            case 'off'
                V='on';
        end
        set(x,'visible',V)
    case 'style'
        s=get(x(1),'linestyle');
        switch lower(s)
            case 'none'
                LineStyle='-';
            otherwise
                LineStyle='none';
        end
        set(x,'linestyle',LineStyle)        
end