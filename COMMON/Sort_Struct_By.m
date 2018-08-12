function P=Sort_Struct_By(P,B)
A=getfield(P,B);
if ~isstr(A{1})
    B=A;A=[];
    for k=1:length(B)
        A(k)=B{k};
    end
end
[X,I]=sort(A);

FN=fieldnames(P);
for k=1:length(FN)    
    p=P.(FN{k});
    p=p(I);
    P.(FN{k})=p;
end