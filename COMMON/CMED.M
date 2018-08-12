function y=cmed(x,s)
flag=0;
[m,l]=size(x);
if m<s
    flag=1;
    x=x';    
    [m,l]=size(x);
end
zs=zeros(fix(s/2),l);
y=emmed1([zs;x;zs],s);
y(y==0)=zeros(size(y(y==0)));
if flag==1
    y=y';
end
clear emmed1
return
[m,l]=size(x);
if m<l
	x=x';
end;
zs=zeros(1,fix(s/2));
	y=mmed1([zs x' zs],s)';
	y(y==0)=zeros(size(y(y==0)));
if m<l
	y=y';
end;
clear mmed1