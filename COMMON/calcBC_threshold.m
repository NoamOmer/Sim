function T=calcBC_threshold(varargin)
% T=calcbc_threshold(Len,Nfft,Olap,Sel,p,win)
v=varargin;
if isstr(v{4})
    [m,Nr0]=size(vec2mat((1:v{1})',v{2},0));
    [m,Nr]=size(vec2mat((1:v{1})',v{2},v{3}));
    v(1:3)=[];
else
    Nr=v{1};Nr0=v{2};v(1:2)=[];
end
sqrN=sqrt(Nr);
logN=log(Nr);
logN2=logN^2;
NN=sqrt(Nr0/(1+(Nr0/Nr)^2));
NN1=sqrt(Nr0/(1/sqrt(2)+(Nr0/Nr)^2));

H.hanning=[1.373*sqrN+0.219,    1.415*sqrN-0.045,    -.016*logN2+.239*logN+2.311,    2.58/sqrN];
H.bartlett=[1.37*sqrN+.23,    1.413*sqrN-0.039,    -.015*logN2+.236*logN+2.319,    2.59/sqrN];
H.kaiser=[1.351*sqrN+0.356,    1.418*sqrN-0.055,    -.015*logN2+.234*logN+2.322,    .92/Nr+2.29/sqrN+0.02];
H.rectangle=[1.982*NN+0.094,    1.83*NN1-0.03,    Nr/(0.314*Nr+0.955),    2.59/sqrN];

KP.hanning=[1.403*sqrN-0.015,    1.415*sqrN-0.044,    -.016*logN2+.238*logN+2.315,    2.58/sqrN];
KP.bartlett=[1.399*sqrN,    1.413*sqrN-0.037,    -.015*logN2+.232*logN+2.328,    2.59/sqrN];
KP.kaiser=[1.392*sqrN+0.053,    1.418*sqrN-0.055,    -.015*logN2+.23*logN+2.333,    2.57/sqrN];
KP.rectangle=[2.012*NN-0.093,    1.83*NN1-0.04,    Nr/(0.314*Nr+0.957),    2.59/sqrN];


thr=lower(v{1});v(1)=[];
p=v{1};v(1)=[];
win=v{1};v(1)=[];
if isempty(v),
    method='haubrich';
else
    method=v{1};
end
if ~isempty(v)
    method=lower(v{1});
end
switch method
    case {'haubrich',1}
        A=H;
    otherwise
        A=KP;
end
%disp(['Nr=',num2str(Nr),'--- Nr0=',num2str(Nr0)])
L=getfield(A,lower(win));
switch thr
    case 'amplitude'
        b=L(1);
        T=sqrt(chi2inv(p,2))/b;
    case 'phase'
        b=L(2);
        T=sqrt(-2*log(1-p))/b;
    case 'vop'
        M=L(3);S=L(4);
        T=(-sqrt(2)*S)*erfcinv(2*p)+M;
        %T=norminv(p,M,S);
        return
end
% sqrt((-2*bk .^ 2) .* log(1 - pk));