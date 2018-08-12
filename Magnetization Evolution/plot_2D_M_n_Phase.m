
function plot_2D_M_n_Phase(M_init,M,M_title_str,Phase_title_str,context)
set_context;

if (~DEBUG_FLAG)  return;  end;

if (~isempty(M_title_str))
	M_vec  = [M{:}];
	Mx_vec = M_vec(1:3:length(M_vec));
	My_vec = M_vec(2:3:length(M_vec));
	Mz_vec = M_vec(3:3:length(M_vec));
	Mx_mat = reshape(Mx_vec,Nx,Ny);
	My_mat = reshape(My_vec,Nx,Ny);
	Mz_mat = reshape(Mz_vec,Nx,Ny);
	Mxy  = sqrt(Mx_mat.^2 + My_mat.^2);
	Mxyz = sqrt(Mx_mat.^2 + My_mat.^2 + Mz_mat.^2);

	% Calculate the magnetization density factor per-region: Z & XY
	M_init_F = sum(sum(M_init)) / sum(sum(M_init));
	Mxy_F    = sum(sum(Mxy))    / sum(sum(M_init));
	Mxyz_F   = sum(sum(Mxyz))   / sum(sum(M_init));
	Mz_mat_F = sum(sum(Mz_mat)) / sum(sum(M_init));

	fh = figure;
	s1 = subplot(2,2,1);
	imagesc(sample_x_axis,sample_y_axis,transpose(M_init));  title(sprintf('Initial M_z(x,y) [F=%2.2f]',M_init_F));
	xlabel('x [cm]');  ylabel('y [cm]');  set(gca,'Ydir','normal'); set(s1,'CLim',[-1.1 1.1]); colorbar;
	s2 = subplot(2,2,2);
	imagesc(sample_x_axis,sample_y_axis,transpose(Mxy));     title(sprintf('M_x_y(x,y)       [F=%2.2f]',Mxy_F   ));
	xlabel('x [cm]');  ylabel('y [cm]');  set(gca,'Ydir','normal'); set(s2,'CLim',[-1.1 1.1]); colorbar;
	s3 = subplot(2,2,3);
	imagesc(sample_x_axis,sample_y_axis,transpose(Mxyz));    title(sprintf('M_{xyz}(x,y)     [F=%2.2f]',Mxyz_F  ));
	xlabel('x [cm]');  ylabel('y [cm]');  set(gca,'Ydir','normal'); set(s3,'CLim',[-1.1 1.1]); colorbar;
	s4 = subplot(2,2,4);
	imagesc(sample_x_axis,sample_y_axis,transpose(Mz_mat));  title(sprintf('M_z(x,y)         [F=%2.2f]',Mz_mat_F));
	xlabel('x [cm]');  ylabel('y [cm]');  set(gca,'Ydir','normal'); set(s4,'CLim',[-1.1 1.1]); colorbar;
% 	set(s2,'CLim',get(s1,'CLim'));
% 	set(s3,'CLim',get(s1,'CLim'));
% 	set(s4,'CLim',get(s1,'CLim'));
    figure(fh);	text(-5,8,M_title_str,'FontWeight','Bold');
end;

if (~isempty(Phase_title_str))
	figure;
	subplot(2,2,1);
	Phi_M2 = unwrap(angle(Mx_mat + 1i*My_mat),[],2); % DIM=2 unwrapping along rows
	imagesc(sample_x_axis,sample_y_axis,transpose(Phi_M2));
	title([Phase_title_str ' (y axis unwrapped)']); xlabel('x-axis [cm]'); ylabel('y-axis [cm]'); set(gca,'Ydir','normal');

	subplot(2,2,2);
	Phi_M1 = unwrap(angle(Mx_mat + 1i*My_mat),[],1); % DIM=1 unwrapping along cols
	imagesc(sample_x_axis,sample_y_axis,transpose(Phi_M1));
	title([Phase_title_str ' (x axis unwrapped)']); xlabel('x-axis [cm]'); ylabel('y-axis [cm]'); set(gca,'Ydir','normal');

    subplot(2,2,3);
    plot(sample_y_axis,transpose(Phi_M2(round(length(Phi_M2(:,1))/2),:)),'.-');
	title([Phase_title_str ' (y axis unwrapped (profile))']); xlabel('y-axis [cm]'); ylabel('Phase [rad]');

    subplot(2,2,4);
    plot(sample_x_axis,transpose(Phi_M1(:,round(length(Phi_M1(1,:))/2))),'.-');
	title([Phase_title_str ' (x axis unwrapped (profile))']); xlabel('x-axis [cm]'); ylabel('Phase [rad]');
end;

return;

