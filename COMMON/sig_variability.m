function y=sig_variability(x,N);
if iseven(N)
    N=N+1;
end
y=comp_sigvariability(x,N);
clear comp_sigvariability