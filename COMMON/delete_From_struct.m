function P=delete_From_struct(P,l);
if isempty(l),return,end
S=fieldnames(P);
for k=1:length(S)    
    p=getfield(P,S{k});
    p(l)=[];
    P=setfield(P,S{k},p);
end
