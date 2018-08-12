function P1=Get_Struct_line(P,lineNo)
P1=P;
S=fieldnames(P);
for k=1:length(S)
    p=getfield(P,S{k});
    P1=setfield(P1,S{k},p(lineNo));
end