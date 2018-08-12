function P=add_struct(P,varargin);
S=fieldnames(P);
if isstruct(varargin{1})
    for l=1:length(varargin)
        for k=1:length(S)    
            p=[];
            p1=getfield(P,S{k});
            if iscell(p1)
                p=p1;
                p{end+1}=getfield(varargin{l},S{k});
            else
                p{1}=p1;
                p{end+1}=getfield(varargin{l},S{k});
            end
            P=setfield(P,S{k},p(:));
        end
    end
else
    if length(varargin)~=length(S),
        errordlg('You must provide a struct variable P + variables for each filed of P');
        return 
    end
    for k=1:length(S)    
        p=getfield(P,S{k});
        p{end+1}=varargin{k};
        P=setfield(P,S{k},p(:));
    end
end