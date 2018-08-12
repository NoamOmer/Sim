function SaveToSpat(X,Xname,Xfile)
if iscell(X),X=X{1};end
if iscell(Xname),Xname=Xname{1};end
if iscell(Xfile),Xfile=Xfile{1};end
fig1=findobj(0,'name','SPAT-Signal Preparation and Analysis Toolbox');
a=get(fig1,'userdata');
ena=1;
K=length(a.varName)+1;
for l=1:K-1
    if strmatch(a.varName{l},Xname)
        b=questdlg('Variable already exists in SPAT environment, replace it?'...
            ,' ','Yes','No','Yes');
        if strcmp(b,'Yes'),ena=1;K=l;else,ena=0;end
    end
end
if ena==1
    a.varName{K}=Xname;
    a.varData{K}=X;
    a.varFile{K}=Xfile;    
end    
set(fig1,'userdata',a);
