function f=isdouble(X)
classX=class(X);
f=0;
if strcmp(lower(classX),'double')
    f=1;
end