
% Prepare structure for 'evolve_M_n_acq1Dgrad_t_CPP' C function
% See the function documentation for parameters description.
function st = set_evolve_M_n_acq1Dgrad_t_CPP_struct(M_init,M0z,B_eff_rot,dt,dB0z,zGrad,z_axis,RH_flag,...
                                                    acq_flag,inhomo_flag,relax_flag,T1,T2,context)

if (nargin ~= 14)
	error('Wrong number of parameters');
end;

st.M_initx     = transpose(M_init   (:,1));    % 1
st.M_inity     = transpose(M_init   (:,2));    % 2
st.M_initz     = transpose(M_init   (:,3));    % 3
st.M0z         = M0z;                          % 4
st.B_eff_rotx  = transpose(B_eff_rot(:,1));    % 5
st.B_eff_roty  = transpose(B_eff_rot(:,2));    % 6
st.B_eff_rotz  = transpose(B_eff_rot(:,3));    % 7
st.dB0z        = dB0z;                         % 8
st.zGrad       = zGrad;                        % 9
st.z_axis      = z_axis;                       % 10
st.dt          = dt;                           % 11
st.RH_flag     = RH_flag;                      % 12
st.acq_flag    = acq_flag;                     % 13
st.inhomo_flag = inhomo_flag;                  % 14
st.relax_flag  = relax_flag;                   % 15
st.T1          = T1;                           % 16
st.T2          = T2;                           % 17
st.context     = context;                      % 18

return;

