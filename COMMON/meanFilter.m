function y=meanFilter(x,L);
if rem(L,2)==0,
    L=L+1;
    warning('Filter length should be odd integer, increasing it by 1');
end
[m,n]=size(x);
if m<L
    x=x';
end
y=comp_meanfilter(x,L);
if m<L
    y=y';
end
clear comp_meanfilter