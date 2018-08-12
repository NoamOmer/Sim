
function M_final = T2_relaxation(M_init, t, T2)

T2_vec  = [exp(-t/T2), exp(-t/T2), 1];

M_final = T2_vec .* M_init;

return;

