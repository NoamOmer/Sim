function M=divide(X,Y)
M=X;
for k=1:length(Y)
    M(k,:)=X(k,:)./Y(k);
end