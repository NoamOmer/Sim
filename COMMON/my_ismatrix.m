function f=my_ismatrix(X)
[m,n]=size(X);
f=1;
if m==1 | n==1
    f=0;
end