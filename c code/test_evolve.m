% ------------------------------
%  Test evolve_M_n_acq1Dgrad_t.c
% ------------------------------
clear all;
mex evolve_M_n_acq1Dgrad_t.c

B_rot = [0,0,0];
dt = 2e-6;
t = 0:2e-6:12e-3;
z = -1.5:2E-3:1.5;
% t = 0:dt:20e-6;
% z = -1.5E-2:2E-3:1.5E-2;

st.M_initz     = [zeros(1,5), ones(1,length(z)-10), zeros(1,5)];
st.M_inity     = zeros(1,length(st.M_initz));
st.M_initx     = zeros(1,length(st.M_initz));
st.B_eff_rotx  = ones(1,length(t)) * B_rot(1);
st.B_eff_roty  = ones(1,length(t)) * B_rot(2);
st.B_eff_rotz  = ones(1,length(t)) * B_rot(3);
st.dB0z        = z*0;
st.zGrad       = ones(1,length(t))*2;
st.z_axis      = z;
st.dt          = dt;
st.RH_flag     = 1;
st.acq_flag    = 1;
st.inhomo_flag = 1;

tic;
[Mx,My,Mz,sig_real,sig_imag] = evolve_M_n_acq1Dgrad_t(st);
tm = toc;
disp(sprintf('Total simulation time: %f [sec]',tm));
disp(sprintf('z len = %f\nt len = %f',length(z),length(t)));

% ----------------
%  Test evolve_M.c
% ----------------
clear all;
mex evolve_M.c

z = -1.5:2E-3:1.5;
% z = -1.5E-2:2E-3:1.5E-2;

st.M_initz     = [zeros(1,5), ones(1,length(z)-10), zeros(1,5)];
st.M_inity     = zeros(1,length(st.M_initz));
st.M_initx     = zeros(1,length(st.M_initz));
st.Tev         = 3e-3;
st.dB0z        = z*0;
st.zGrad       = 2;
st.z_axis      = z;
st.RH_flag     = 1;
st.inhomo_flag = 1;

tic;
[Mx,My,Mz] = evolve_M(st);
tm = toc;
disp(sprintf('Total simulation time: %f [sec]',tm));
disp(sprintf('z len = %f\n',length(z)));
