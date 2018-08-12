function RRsignal=comp_rr(Tr,RR,rrSF)
if nargin==1
    rrSF=10;
    RR=diff(Tr);RR=[RR(1);RR];
elseif nargin==2
    if size(Tr)==size(RR)
        rrSF=10;
    else
        rrSF=RR;
        RR=diff(Tr);RR=[RR(1);RR];
    end           
end
if (Tr(2)-Tr(1))<RR(1)
    SF=RR(1)/(Tr(2)-Tr(1))
    RR=RR./SF;
end    
RRsignal=comp_rr(Tr,RR,rrSF);
clear comp_rr