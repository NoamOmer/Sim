function HR=make_hr(Tr,RR,hrSF)
if nargin==1
    hrSF=10;
    RR=diff(Tr);
    RR=[RR(1);RR];
elseif nargin==2
    if length(Tr(:))==length(RR(:))
        hrSF=10;
    else
        hrSF=RR;
        RR=diff(Tr);
        RR=[RR(1);RR];
    end           
end
% if (Tr(2)-Tr(1))<RR(1)
%     SF=RR(1)/(Tr(2)-Tr(1))
%     RR=RR./SF;
% end    
clear comp_hr
HR=comp_hr(Tr,RR,hrSF);
clear comp_hr