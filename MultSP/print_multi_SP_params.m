
function print_multi_SP_params(GePE, Tp, Ly, NPE, Nsp)

set_globals;
nn = 0:5;
OPnp = round((gammaHz * GePE * Tp * Ly) ./ (2*nn + 1));
declare_stage('Multi-SP paraemters optimization');
disp(sprintf('Optimum    np  :  %3.0f  %3.0f  %3.0f  %3.0f ...',OPnp(1),OPnp(2),OPnp(3),OPnp(4)));
disp(sprintf('Simulation np  :  %3.0f                         ',2*NPE                          ));

nMy = 2*NPE;
dy  = Ly / (Nsp*nMy);
Dk  = 1/dy;
k0_tag = gammaHz*GePE*Tp;
k0_tag = mod(k0_tag,Dk);
disp(sprintf('\n  Dk = [%3.1f ... %3.1f]\n  k0_tag = %3.1f (filetered-out echo location)\n',Dk/2,+Dk/2,k0_tag));

return;
