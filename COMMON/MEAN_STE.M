function varargout=mean_ste(X)
if nargout<2
    m=compMeanStd(X);
    varargout{1}=[m(1,:);m(2,:)./sqrt(size(X,1))];
else
    [me,st]=compMeanStd(X);
    varargout{1}=me;
    varargout{2}=st./sqrt(size(X,1));
end
clear compMeanStd