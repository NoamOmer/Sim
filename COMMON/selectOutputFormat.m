function F=selectOutputFormat;
a=listdialog('listsize',[150,80],'liststring',[{'Mat'};{'Acq'};{'Txt'};{'int8'};{'int16'};{'float'};{'double'}]...
    ,'selectionmode','single','name','FORMAT','promptstring','Select output FORMAT');
if a==1, 
    F='mat';
elseif a==2, 
    F='acq';
elseif a==3, 
    F='txt';
elseif a==4, 
    F='int8';
elseif a==5, 
    F='int16';
elseif a==6, 
    F='float';
elseif a==7, 
    F='double';
end