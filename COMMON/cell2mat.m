function A=cell2mat(P)
A=[];
for k=1:size(P,2)
    B=[];
    for l=1:size(P,1)
        B=[B;P{l,k}];
    end
    A=[A,B];
end