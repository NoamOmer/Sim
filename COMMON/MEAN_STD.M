function varargout=mean_std(X)
if nargout<2
    m=compMeanStd(X);
    varargout{1}=m;
else
    [me,st]=compMeanStd(X);
    varargout{1}=me;
    varargout{2}=st;
end
clear compMeanStd