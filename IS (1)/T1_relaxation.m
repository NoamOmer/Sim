% T1 relaxation function.
% Note that M0z does not denote the initial magnetization but rather 
% the magnetization at equiibrium (pointing along the field axis).

function   M_final  = T1_relaxation(M_init, M0z, t, T1)

M_final = [M_init(1), M_init(2), M_init(3)*exp(-t/T1) + M0z*(1 - exp(-t/T1))];

return;

