function X=vec2mat(x,N,lap);
if nargin<3
    h=errordlg('vec2mat requires 3 inputs','Vec2mat Error','modal');
    waitfor(h);
    return
end
if iscell(x)
    for k=1:length(x)
        X{k}=performFCN(x{k},N,lap);
    end
else
    X=performFCN(x,N,lap);
end

function X=performFCN(x,N,lap);
if N>length(x)
    h=errordlg('Matrix raws cannot exceed vector length','Vec2mat Error','modal');
    waitfor(h);
    return
end
if lap>N
    h=errordlg('Overlap should be smaller than matrix length','Vec2mat Error','modal');
    waitfor(h);
    return
end
X=vector2matrix(x,N,lap);
clear vector2matrix