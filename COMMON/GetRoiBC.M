function K=GetRoiBC(N,k)
if nargin==1    
    K=compRoi(N);
else
    K=compRoi(N,k);
end
clear compRoi