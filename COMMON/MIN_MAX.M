function varargout=min_max(X,opt)
if nargin==1,opt='removeNaN';end
switch opt
case 'removeNaN'
    mi=min(X);
    ma=max(X);
otherwise
    [mi,ma]=compMinMax(X);
    clear compMinMax
end

if nargout<2
    varargout{1}=[mi;ma];
else
    varargout{1}=mi;
    varargout{2}=ma;
end