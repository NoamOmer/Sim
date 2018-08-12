function generate_signals
t=[1000,1600,2600;
    1100,1800,3000;
    1200,1800,2600;
    1300,1900,3000;
    1400,2100,2900;
    1500,2100,3000];
val=[50,70,60;
    50,70,50;
    60,70,60;
    80,60,80;
    70,40,70;
    90,70,90];
for k=1:6
    LF=generate(t(k,:),val(k,:));
    save(['c:\data\sim\d',num2str(k),'\HR_integrals'],'LF')
end

function x=generate(t,val)
x(1:t(1))=val(1);
T=t(1):t(2);
tau=(t(2)-t(1))/10;
x(T)=val(1)+(val(2)-val(1)).*(1-exp(-(T-t(1))/tau));
T=t(2):t(3);
tau=(t(3)-t(2))/10;
x(T)=(val(2)-val(3)).*exp(-(T-t(2))/tau)+val(3);
x=x+randn(size(x));