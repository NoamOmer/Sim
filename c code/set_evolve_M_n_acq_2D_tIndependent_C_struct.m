
% Prepare structure for 'evolve_M_n_acq_2D_tIndependent_C' C function
% See the function documentation for parameters description.
function st = set_evolve_M_n_acq_2D_tIndependent_C_struct(M_init,M0z,B_eff_rot,dB0_XY,Grad,Tevolution,x_axis,y_axis,...
                                                          inhomo_flag,relax_flag,plot_flag,T1,T2,context)
if (nargin ~= 14)
	error('Wrong number of parameters');
end;

	
if (iscell(M_init))
	% Convert cell matrix 1D vector: {1,1},{2,1},{3,1}...{Nx,1},{1,2},{2,2}...{Nx,2},...{Nx,Ny}
	% where each of the above cells is a triplet of Mx,My,Mz values
	M_vec  = [M_init{:}];

	% Convert 1D vector into x,y,x vectors: 
	Mx_vec = M_vec(1:3:length(M_vec));  % Mx values of (x,y) locations: {1,1},{2,1},{3,1}...{Nx,1},{1,2},{2,2}...{Nx,2},...{Nx,Ny}
	My_vec = M_vec(2:3:length(M_vec));
	Mz_vec = M_vec(3:3:length(M_vec));
else
	Mx_vec = M_init(1,:);
	My_vec = M_init(2,:);
	Mz_vec = M_init(3,:);
end;

dB0_XY = transpose(dB0_XY(:));

st.M_initx     = Mx_vec;                       % 1   vec     2D XY matrix converted to 1D vector
st.M_inity     = My_vec;                       % 2   vec     2D XY matrix converted to 1D vector
st.M_initz     = Mz_vec;                       % 3   vec     2D XY matrix converted to 1D vector
st.M0z         = M0z;                          % 4
st.B_eff_rot   = B_eff_rot;                    % 5   vec
st.dB0_XY      = dB0_XY;                       % 6   vec     2D XY matrix converted to 1D vector
st.Ge          = Grad;                         % 7   vec
st.x_axis      = x_axis;                       % 8   vec
st.y_axis      = y_axis;                       % 9   vec
st.Tevolution  = Tevolution;                   % 10   scalar
st.inhomo_flag = inhomo_flag;                  % 11  scalar
st.relax_flag  = relax_flag;                   % 12  scalar
st.plot_flag   = plot_flag;                    % 13  scalar
st.T1          = T1;                           % 14  scalar
st.T2          = T2;                           % 15  scalar
st.context     = context;                      % 16  string

return;

