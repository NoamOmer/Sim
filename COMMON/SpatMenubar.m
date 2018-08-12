function SpatMenubar(action,fig)
if nargin<2
    fig=gcbf;
end
if ~strcmp(get(fig,'type'),'figure')
    fig=gcbf;
end  
% MainBar=findHiddenObjects(fig,'tag','FigureToolBar');
% if isempty(MainBar)
%     action='Create';
% else
%     action='Refresh';
% end    
switch action
    case 'Refresh'
        MainBar=findHiddenObjects(fig,'userdata','FigureToolBar');
        a=get(MainBar,'children');
        for k=1:length(a)
            if strcmp(lower(get(a(k),'state')),'on')
                set(a(k),'state','off')
            end
        end
        SubMenu=findHiddenObjects(fig,'tag','subMenu');
        set(SubMenu,'visible','off')
        return
    case 'Create'
        if strcmp(lower(get(gcf,'menubar')),'none')
            return
        end
        temp=uicontrol('parent',fig); %Create temporary uicontrol
        domymenu('menubar','toggletoolbar',fig) %activate the toolbar
        delete(temp)% delete the uicontrol
        set(fig,'menu','figure')
        deleteHiddenObjects(findHiddenObjects(fig,'type','uitoolbar','Tag','subMenu'))
        MainBar=findHiddenObjects(fig,'type','uitoolbar');
        set(MainBar,'userdata','FigureToolBar')
        a=get(MainBar,'children');%get(a,'CreateFcn')
%         set(a,'CreateFcn',';')
        subMenu=copyobj(MainBar,fig);
        set(subMenu,'visible','off','tag','subMenu');
        deleteHiddenObjects(get(subMenu,'children'))

        b=findHiddenObjects(MainBar,'type','uitoggletool');
        for k=1:length(b)
            S=[get(b(k),'CreateFcn'),get(b(k),'tag'),get(b(k),'tooltip')];
            if ~isempty(findstr(lower(S),'zoom'))
                delete(b(k))
            else
                %set(b(k),'clickedcallback','','oncallback','','offcallback','','CreateFcn','')
            end
        end
        b=b(1);
        zoomBtn=copyobj(b,MainBar);
        autoscaleBtn=copyobj(zoomBtn,MainBar);
        show_pointer=copyobj(zoomBtn,MainBar);
        add_line=copyobj(zoomBtn,MainBar);
        shapes=copyobj(zoomBtn,MainBar);
%         get(zoomBtn)
        C=get_cdata('ShowLoc');
        set(show_pointer,'clickedcallback','','cdata',C,'tooltipstring','Track pointer','oncallback',...
            'set(gcf,''windowbuttonmotion'',''show_pointer_location(1);'')','tag','show_pointer',...
            'offcallback','set(gcf,''windowbuttonmotion'',''show_pointer_location(0);'')')
%******************************************************************************************
% Zoom
%******************************************************************************************
        f=[];
        a=zoomBtn;
        C=get_cdata('zoominXY');    
        f(1)=copyobj(a,subMenu);
        set(f(1),'tooltipstring','Zoom current axis','tag','zoomGcaXY','state','on','cdata',C,...
            'oncallback','','offcallback','')
        
        C=get_cdata('zoominX');    
        f(2)=copyobj(f(1),subMenu);
        set(f(2),'tooltipstring','Zoom current axis','tag','zoomGcaX','separator','off','cdata',C,...
            'state','off')
        
        C=get_cdata('zoominY');    
        f(3)=copyobj(f(1),subMenu);
        set(f(3),'tooltipstring','Zoom current axis','tag','zoomGcaY','separator','off','cdata',C,...
            'state','off')
        
        C=get_cdata('zoominXYX');        
        f(4)=copyobj(f(1),subMenu);
        set(f(4),'tooltipstring','Zoom current axis','tag','zoomXYX','separator','off','cdata',C,...
            'state','off')
        
        C=get_cdata('zoominXX');    
        f(5)=copyobj(f(1),subMenu);
        set(f(5),'tooltipstring','Zoom current axis','tag','zoomXX','separator','off','cdata',C,...
            'state','off')
        
        set(f,'userdata',zoomBtn,'clickedcallback','SpatMenubar(''toggleSubmenuObject'',gcf)')
        
        C=get_cdata('zoom');
        set(zoomBtn,'cdata',C,'clickedcallback','','tooltipstring','Zoom In All Axes',...
            'tag','xZoomInAll','oncallback','SpatMenubar(''ZoomMenuOn'',gcbf)',...
            'offcallback','SpatMenubar(''ZoomMenuOff'',gcbf)','separator','On')
        
%******************************************************************************************
 % Auto Scale submenu
%******************************************************************************************               
       a=findHiddenObjects(fig,'type','uipushtool');a=a(1);
       C=get_cdata('autoscaleGcaXY');    
       f=copyobj(a,subMenu);
       set(f,'separator','on','clickedcallback','AutoScale(''GcaXY'',gcf)','tooltipstring',...
            'AutoScale current axis','cdata',C,'userdata',autoscaleBtn)
        
        C=get_cdata('autoscaleGcaX');    
        f=copyobj(a,subMenu);
        set(f,'clickedcallback','AutoScale(''GcaX'',gcf)','tooltipstring',...
            'AutoScale X Gca','cdata',C,'userdata',autoscaleBtn)
        
        C=get_cdata('autoscaleGcaY');    
        f=copyobj(a,subMenu);
        set(f,'clickedcallback','AutoScale(''GcaY'',gcf)','tooltipstring',...
            'AutoScale Y Gca','cdata',C,'userdata',autoscaleBtn)
        
        C=get_cdata('autoscaleXY');
        f=copyobj(a,subMenu);
        set(f,'clickedcallback','AutoScale(''XY'',gcf)','tooltipstring',...
            'AutoScale All axes','cdata',C,'userdata',autoscaleBtn)
        
        C=get_cdata('autoscaleXX');    
        f=copyobj(a,subMenu);
        set(f,'clickedcallback','AutoScale(''XX'',gcf)','tooltipstring',...
            'AutoScale X All Axes','cdata',C,'userdata',autoscaleBtn)
        
        C=get_cdata('autoscaleXYX');    
        f=copyobj(a,subMenu);
        set(f,'clickedcallback','AutoScale(''XYX'',gcf)','tooltipstring',...
            'AutoScale X All & Y Gca','cdata',C,'userdata',autoscaleBtn)
        
        C=get_cdata('autoscaleYY');    
        f=copyobj(a,subMenu);
        set(f,'clickedcallback','AutoScale(''YY'',gcf)','tooltipstring',...
            'AutoScale Y All axes','cdata',C,'userdata',autoscaleBtn)
        
        C=get_cdata('autoscaleGcaXY');    
        set(autoscaleBtn,'clickedcallback','','tooltipstring','Auto Scale','cdata',C,...
            'tag','AutoScale','oncallback','SpatMenubar(''AutoScaleMenuOn'',gcbf)',...
            'offcallback','SpatMenubar(''AutoScaleMenuOff'',gcbf)','separator','off')
%******************************************************************************************
% Add lines
%******************************************************************************************        
        g=[];
        C=get_cdata('AddXline');    
        g(1)=copyobj(f,subMenu);
        set(g(1),'tooltipstring','Add horizontal line (gca)','tag','gcaHoriz','separator','on','cdata',C)
        
        C=get_cdata('AddXXline');    
        g(2)=copyobj(f,subMenu);
        set(g(2),'tooltipstring','Add horizontal line','tag','Horiz','separator','on','cdata',C)
        
        C=get_cdata('AddYline');    
        g(3)=copyobj(f,subMenu);
        set(g(3),'tooltipstring','Add vertical line','tag','Vert','separator','on','cdata',C)
        
        set(g,'clickedcallback','toolbar_AddLine','userdata',add_line)
        
        C=get_cdata('Addline');    
        set(add_line,'clickedcallback','','tooltipstring','Add lines','cdata',C,...
            'tag','AddLine','oncallback','SpatMenubar(''AddLineOn'',gcbf)',...
            'offcallback','SpatMenubar(''AddLineOff'',gcbf)','separator','off')

%******************************************************************************************
% Add auto shapes
%******************************************************************************************
        g=[];        
        C=get_cdata('rect');    
        g(1)=copyobj(f,subMenu);
        set(g(1),'tooltipstring','Add rectangle','tag','rectangle','separator','on','cdata',C)
        
        C=get_cdata('ellipse');    
        g(2)=copyobj(f,subMenu);
        set(g(2),'tooltipstring','Add ellipse','tag','ellipse','separator','off','cdata',C)
        
        C=get_cdata('triangle');    
        g(3)=copyobj(f,subMenu);
        set(g(3),'tooltipstring','Add triangle','tag','triangle','separator','off','cdata',C)
        
        set(g,'clickedcallback','Addshapes','userdata',shapes)
        
        C=get_cdata('Shapes');    
        set(shapes,'clickedcallback','','tooltipstring','Autoshapes','cdata',C,...
            'tag','Shapes','oncallback','SpatMenubar(''ShapesOn'',gcbf)',...
            'offcallback','SpatMenubar(''ShapesOff'',gcbf)','separator','off')
        
        a=[zoomBtn,autoscaleBtn,add_line,shapes];
        set(a,'userdata',a);
        return    
    case 'ZoomMenuOn'
        setFamilyOff(gcbo);
        turnChildrenOn(gcbo);
        zoomallaxes('zoomin',gcf)
    case 'ZoomMenuOff'
        xzoomallaxes(gcf,'off')
    case 'toggleSubmenuObject'
        toggleFamily(gcbo);
    case {'AutoScaleMenuOn','AddLineOn','ShapesOn'}
        setFamilyOff(gcbo);
        turnChildrenOn(gcbo);
    case 'AutoScaleMenuOff'
    case 'AddLineOff'   
    case 'ShapesOff'
end

SubMenu=findHiddenObjects(fig,'tag','subMenu');
MainBar=findHiddenObjects(fig,'userdata','FigureToolBar');
if get(gcbo,'parent')~=SubMenu
    h=get(gcbo,'userdata');
    l=findHiddenObjects(h,'state','On');
    if isempty(l),
        set(SubMenu,'visible','Off'),
    else,  
        set(SubMenu,'visible','On');
    end
end



function setFamilyOff(h)
a=get(h,'userdata');
l=find(a==h);
a(l)=[];
set(a,'state','Off')

function turnChildrenOn(h)
a=findHiddenObjects(gcbf,'tag','subMenu');
% set(a,'visible','on')
b=findHiddenObjects(a,'parent',a);
set(b,'enable','off')
l=findHiddenObjects(b,'userdata',h);
set(l,'enable','on')

function  toggleFamily(h)
a=get(h,'parent');
u=get(h,'userdata');
l=findHiddenObjects(a,'userdata',u);
set(l,'state','off')
set(h,'state','on')


function C=get_cdata(opt);
N=16;
R=NaN*ones(N);
G=NaN*ones(N);
B=NaN*ones(N);
switch opt
    case 'zoom'
        for k=1:5,  R(17-k,k)=0; R(16-k,k)=0;        end
        R(5:8,4)=0;R(3:4,5)=0;R(2,6:7)=0;R(1,8:11)=0;R(2,12:13)=0;R(3:4,14)=0;
        R(5:8,15)=0;R(9:10,5)=0;R(11,6:7)=0;R(12,8:11)=0;R(11,12:13)=0;R(9:10,14)=0;
        R(4:9,9:10)=1;R(6:7,7:12)=1;
    case 'zoominXY'
        R(2:15,2:3)=0;
        R(14:15,2:15)=0;
        % X
        R(4:5,5)=0;    R(6:7,6)=0;    R(8:9,7)=0;   R(10:11,8)=0;   R(12:13,9)=0;    
        R(4:5,9)=0;    R(6:7,8)=0;    R(8:9,7)=0;   R(10:11,6)=0;   R(12:13,5)=0;    
        % Y
        R(4:5,11)=0;   R(6:7,12)=0;   R(8:9,13)=0;  R(10:11,12)=0;  R(12:13,11)=0;
        R(4:5,15)=0;   R(6:7,14)=0;
        G=R;B=R;
    case 'zoominX'
        R(2:15,2:3)=0;
        R(14:15,2:15)=0;
        % Y
        R(4:5,7)=0;    R(6:7,8)=0;    R(8:9,9)=0;   R(10:11,10)=0;   R(12:13,11)=0;
        R(4:5,11)=0;   R(6:7,10)=0;   R(8:9,9)=0;   R(10:11,8)=0;    R(12:13,7)=0;    
        G=R;B=R;
    case 'zoominY'
        R(2:15,2:3)=0;
        R(14:15,2:15)=0;
        % X
        R(4:5,7)=0;    R(6:7,8)=0;    R(8:9,9)=0;       
        R(4:5,11)=0;   R(6:7,10)=0;   R(8:9,9)=0;   R(10:11,8)=0;    R(12:13,7)=0;    
        G=R;B=R;
    case 'zoominXX'
        % Top axes
        R(2:7,2)=0;
        R(7,2:15)=0;
        R(2,7)=0;    R(3,8)=0;   R(4,9)=0;    R(5,10)=0;    R(6,11)=0;    
        R(2,11)=0;   R(3,10)=0;  R(4,9)=0;    R(5,8)=0;     R(6,7)=0;    
        % Bottom Axes
        R(9:15,2)=0;
        R(15,2:15)=0;
        R(9,7)=0;    R(10,8)=0;   R(11,9)=0;    R(12,10)=0;    R(13,11)=0;    
        R(9,11)=0;   R(10,10)=0;  R(11,9)=0;    R(12,8)=0;     R(13,7)=0;    
        G=R;B=R;
    case 'zoominXYX'
        % Top axes
        R(2:7,2)=0;
        R(7,2:15)=0;
        % X
        R(2,4)=0;    R(3,5)=0;   R(4,6)=0;    R(5,7)=0;    R(6,8)=0;    
        R(2,8)=0;    R(3,7)=0;   R(4,6)=0;    R(5,5)=0;    R(6,4)=0;    
        % Y
        R(2,11)=0;    R(3,12)=0;   R(4,13)=0;    
        R(2,15)=0;    R(3,14)=0;   R(4,13)=0;    R(5,12)=0;    R(6,11)=0;    
        % Bottom Axes
        R(9:15,2)=0;
        R(15,2:15)=0;
        R(9,7)=0;    R(10,8)=0;   R(11,9)=0;    R(12,10)=0;    R(13,11)=0;    
        R(9,11)=0;   R(10,10)=0;  R(11,9)=0;    R(12,8)=0;     R(13,7)=0;    
        G=R;B=R;
    case 'autoscaleGcaXY'
        R(2:13,3)=0; R(3,[2,4])=0;R(12,[2,4])=0;
        R(14,4:15)=0;R([13,15],5)=0;R([13,15],14)=0;
        G=R;B=R;
    case 'autoscaleGcaX'
        R(14,4:15)=0;R([13,15],5)=0;R([13,15],14)=0;
        G=R;B=R;
    case 'autoscaleGcaY'
        R(2:13,3)=0; R(3,[2,4])=0;R(12,[2,4])=0;
        G=R;B=R;
    case 'autoscaleXY'
        R(2:6,3)=0; R(3,[2,4])=0;R(5,[2,4])=0;
        R(7,4:15)=0;R([6,8],5)=0;R([6,8],14)=0;
        R(9:13,3)=0; R(10,[2,4])=0;R(12,[2,4])=0;
        R(14,4:15)=0;R([13,15],5)=0;R([13,15],14)=0;
        G=R;B=R;
    case 'autoscaleXYX'
        R(2:6,3)=0; R(3,[2,4])=0;R(5,[2,4])=0;
        R(7,4:15)=0;R([6,8],5)=0;R([6,8],14)=0;
        R(9:13,3)=0; R(10,[2,4])=0;R(12,[2,4])=0;
        G=R;B=R;
    case 'autoscaleYY'
        R(2:6,3)=0; R(3,[2,4])=0;R(5,[2,4])=0;
        R(9:13,3)=0; R(10,[2,4])=0;R(12,[2,4])=0;
        G=R;B=R;
    case 'autoscaleXX'
        R(7,4:15)=0;R([6,8],5)=0;R([6,8],14)=0;
        R(14,4:15)=0;R([13,15],5)=0;R([13,15],14)=0;
        G=R;B=R;
    case 'ShowLoc'
        R(1:3,1)=0;R(1,1:3)=0;
        for k=1:7,      
            R(k,k)=0;
        end
        R(7:end,9:end)=0.7;
        G=R;
        G(8,10)=0;G(9,11)=0;G(10,12)=0;
        G(10,10)=0;G(9,11)=0;G(8,12)=0;
        
        G(12,10)=0;G(13,11)=0;G(14,12)=0;
        G(12,14)=0;G(13,13)=0;G(14,12)=0;G(15,11)=0;G(16,10)=0;
        B=G;
    case 'Addline'
        R(1:7,1)=0;
        R(7,1:16)=0;
        R(9:16,1)=0;
        R(16,1:16)=0;
        G=R;B=R;
        k=[2:2:7];
        R(k,8)=1;
        k=[2:2:16];
        R(13,k)=1;
    case 'AddXline'
        R(1:16,1)=0;
        R(16,1:16)=0;
        G=R;B=R;
        k=2:2:16;
        R(k,8)=1;
    case 'AddXXline'
        R(1:7,1)=0;
        R(7,1:16)=0;
        R(9:16,1)=0;
        R(16,1:16)=0;
        G=R;B=R;
        k=[2:2:6,10:2:15];
        R(k,8)=1;
    case 'AddYline'
        R(1:16,1)=0;
        R(16,1:16)=0;
        G=R;B=R;
        k=2:2:16;
        R(8,k)=1;
    case 'Shapes'
        X=vec2mat(1:N^2,N,0);
        X1=X';
        E=[X(1:6,1);X(1:6,6);X1(1:6,1);X1(1:6,6)];
        R(E)=0;G(E)=0;B(E)=0;E=X(2:5,2:5);R(E)=1;G(E)=0;B(E)=0;
        E=[X1(4:10,10);X(9,5);X(8,6);X(7,7);X(8,8);X(9,9)];
        R(E)=0;G(E)=0;B(E)=0;
        E=[X(9,6:8)';X(8,7)];R(E)=0;G(E)=0;B(E)=1;
        E=[X(10:11,10);X(9,11);X(8,12:13)';X(9,14);X(10:11,15);X(12,14);X(13,12:13)';X(12,11)];
        R(E)=0;G(E)=0;B(E)=0;
        E=[X(10:11,11:14);X(9:12,12:13)'];E=E(:);R(E)=0;G(E)=1;B(E)=1;        
    case 'rect'
        X=vec2mat(1:N^2,N,0);
        X1=X';
        E=[X(2:15,2);X(2:15,15);X1(2:15,2);X1(2:15,15)];
        R(E)=0;G(E)=0;B(E)=0;E=X(3:14,3:14);R(E)=1;G(E)=0;B(E)=0;
    case 'ellipse'
        X=vec2mat(1:N^2,N,0);
        a=7;b=7;
        x=-a:a;
        y=ceil(b/a*sqrt(a*a-x.*x));
        x=[x(:);x(:)];
        y=[y(:);-y(:)];
        X([x+8,y+8])=1;
        for k=1:length(x)
            R(y(k)+8,x(k)+8)=1;
        end
end
C(:,:,1)=R;
C(:,:,2)=G;
C(:,:,3)=B;
