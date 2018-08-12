
% Prepare structure for 'evolve_M_CPP' C function
% See the function documentation for parameters description.
function  st = set_evolve_M_CPP_struct(M,M0z,B_eff_rot,Exc_Delay,dB0z,zGrad,z_axis,RH_flag,...
                                       inhomo_flag,relax_flag,T1,T2,context)
if (nargin ~= 13)
	error('Wrong number of parameters');
end;

st.M_initx     = transpose(M(:,1));         % 1
st.M_inity     = transpose(M(:,2));         % 2
st.M_initz     = transpose(M(:,3));         % 3
st.M0z         = M0z;                       % 4
st.B_eff_rotx  = B_eff_rot(1);              % 5
st.B_eff_roty  = B_eff_rot(2);              % 6
st.B_eff_rotz  = B_eff_rot(3);              % 7
st.Tev         = Exc_Delay;                 % 8
st.dB0z        = dB0z;                      % 9
st.zGrad       = zGrad;                     % 10
st.z_axis      = z_axis;                    % 11
st.RH_flag     = RH_flag;                   % 12
st.inhomo_flag = inhomo_flag;               % 13
st.relax_flag  = relax_flag;                % 14
st.T1          = T1;                        % 15
st.T2          = T2;                        % 16
st.context     = context;                   % 17

return;


