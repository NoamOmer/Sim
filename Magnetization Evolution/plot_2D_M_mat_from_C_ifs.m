
function plot_2D_M_mat_from_C_ifs(Mx,My,Mz,Nx,Ny,n_jumps,context)
set_context;
global plot_2D_M_mat_from_C_ifs_idx;

if isempty(plot_2D_M_mat_from_C_ifs_idx)
	plot_2D_M_mat_from_C_ifs_idx = 1;
else
	plot_2D_M_mat_from_C_ifs_idx = plot_2D_M_mat_from_C_ifs_idx + 1;
end;

if (~DEBUG_FLAG)  return;  end;

% dbstop in plot_2D_M_mat_from_C_ifs at 15
str = sprintf('Calling plot_2D_M_mat_from_C_ifs for the %2.0f time out of %2.0f\n',...
               plot_2D_M_mat_from_C_ifs_idx, n_jumps+1);
disp(str);

M = [transpose(Mx),transpose(My),transpose(Mz)];
M = mat2cell(M,ones(1,length(M)),3);
M = reshape(M,Nx,Ny);

str = sprintf('Debug point %2.0f (%2.0f)\n', plot_2D_M_mat_from_C_ifs_idx, n_jumps+1);
plot_2D_M_n_Phase(M0_sample,M,str,str,context);

% Reset the index if we are finished with this series
if (plot_2D_M_mat_from_C_ifs_idx == n_jumps+1)
	plot_2D_M_mat_from_C_ifs_idx = 0;
end;

pause(0.25);

return;

