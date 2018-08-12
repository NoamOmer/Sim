function K=get_roi(N,k)
if nargin==1    
    K=compRoiBC(N);
else
    K=compRoiBC(N,k);
end
clear compRoiBC