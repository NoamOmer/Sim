function [selection,value] = listdialog(varargin)
v=varargin;
FigureUnits=get(0,'defaultfigureunits');
set(0,'defaultfigureunits','pixels');
[FigName,Prompt,SelectionMode,ListSize,ListString]=get_variables(v);
[selection,value] = listdlg('name',FigName,'promptstring',Prompt,...
    'ListSize',ListSize,'ListString',ListString);
set(0,'defaultfigureunits',FigureUnits);

function [FigName,PromptString,SelectionMode,ListSize,ListString]=get_variables(v);
FigName = '';
SelectionMode = 'multiple';
PromptString ='';
ListString ='';
ListSize = [160 300];
for i=1:2:length(v)
    switch lower(v{i})
        case 'name'
            FigName = v{i+1};
        case 'promptstring'
            PromptString = v{i+1};
        case 'selectionmode'
            SelectionMode=v{i+1};
        case 'listsize'
            ListSize = v{i+1};
        case 'liststring'
            ListString = v{i+1};
    end
end
