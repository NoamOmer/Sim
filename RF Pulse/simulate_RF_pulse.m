clear all; %close all;
tic
context = 'set_simulate_RF_pulse_globals';
set_context;

% ------------------
% Init parameters
% ------------------
% RF file
RootDir = 'D:\PhD\06 Experiments\2DHybrid 2D pulses\Report 3 - mouse in muerto\Seq_params\';
RFname  = 'nbe_chrp180_11809_153617.RF';

% RootDir = 'D:\PhD\Matlab\Simulations\Ver_4\RF Pulse\';
% RFname  = 'test_RF2.RF';
rfwdth  = 4e-3;

% RF power
RF_max_amp_Hz = 5822.9; %1.0522E+001 * 1e3; %5822.9; % [Hz]
RF_max_amp_T  = 1e-4*RF_max_amp_Hz / gammaHz;    % [T]  = 1e-4 * [Hz] / [Hz/Gauss]

% Grad
Ge = [0,5,0];

% Axes
Yie     = -1.5;
Yfe     = +1.5;
Ly      = abs(Yfe-Yie);
Lx      = 2.0;
dx      = 0.05;
dy      = 0.005 * sign(Yfe-Yie);
x_axis  = (-Lx/2): dx : (+Lx/2);
y_axis  = (Yie  ): dy : (Yfe  );
Nx      = length(x_axis);
Ny      = length(y_axis);

% Inhomogeneity
inhomo_flag = 0;
dB0_XY = zeros(Nx,Ny);

% ------------------
% create a sample
% ------------------
M0 = [0,1,0];

M0_sample = 0;
wx   = 0.8;      % square width, relative to image width (normalized to 0..1)
wy   = 0.6;      % square hight, relative to image hight (normalized to 0..1)
smooth_Fx = 5;
smooth_Fy = 5;

sftx = 0.0;      % horizontal shift (normalized to 0..1)
sfty = 0.0;     % vertical   shift (normalized to 0..1)
lx1 = round(((1-wx)/2 - sftx)*Nx);
ly1 = round(((1-wy)/2 - sfty)*Ny);
lx2 = round(wx*Nx);
ly2 = round(wy*Ny);
lx3 = Nx - lx1 - lx2;
ly3 = Ny - ly1 - ly2;
sx = [ones(1,lx1)*1E-4  ,  ones(1,lx2)  ,  ones(1,lx3)*1E-4];
sy = [ones(1,ly1)*1E-4  ,  ones(1,ly2)  ,  ones(1,ly3)*1E-4];
sx = transpose(smooth(sx,smooth_Fx));
sy = transpose(smooth(sy,smooth_Fy));
M0_sample = M0_sample + (transpose(sx))*sy;

M0_sample = M0_sample ./ max(max(M0_sample));
M0_sample = M0_sample + 0.01;

M_init = {};
for x_idx = 1:Nx
	for y_idx = 1:Ny
		M_init{x_idx,y_idx} = M0_sample(x_idx,y_idx) * M0;
	end;
end;

if (DEBUG_FLAG >= 3)
figure;
imagesc(x_axis,y_axis,transpose(M0_sample)); set(gca,'Ydir','normal'); %axis image;
title(sprintf('Initial M_z(x,y)')); xlabel('x [cm]'); ylabel('y [cm]');
set(get(gca,'title') ,'FontWeight','BOLD'); % colorbar;
end;

% -----------------------------------------------
% Load RF
% -----------------------------------------------
rf = load([RootDir RFname]);
rf_amp_norm = (rf(:,2) ./ max(rf(:,2)));
rf = (rf_amp_norm*RF_max_amp_T) .* exp(i*2*pi*rf(:,1)/360);
rf = transpose(rf);
t  = linspace(0,rfwdth,length(rf));

if (DEBUG_FLAG >= 3)
figure;
plot(t*1e3,rf_amp_norm*RF_max_amp_Hz,'.-'); title('Abs RF'  ); xlabel('Time [ms]'); ylabel('Power [Hz]');
figure;
subplot(2,1,1); plot(t*1e3,phase(rf),'.-'); title('Phsae RF (unwrapped)'); xlabel('Time [ms]'); ylabel('Phase');
subplot(2,1,2); plot(t*1e3,angle(rf),'.-'); title('Phsae RF (wrapped)');   xlabel('Time [ms]'); ylabel('Phase');
end;

% -----------------------------------------------
% Apply RF
% -----------------------------------------------
rfx = real(rf);
rfy = imag(rf);
rfz = zeros(1,length(rf));
B_eff = [transpose(rfx),transpose(rfy),transpose(rfz)];
dte   = rfwdth / length(rfx);

% EXCITE
acq_flag  = 0;
plot_flag = 0;

ParamsSt = set_evolve_M_n_acq_2D_C_struct(M_init,B_eff,dB0_XY,Ge,dte,x_axis,y_axis,...
									acq_flag,inhomo_flag,Relax_Flag,plot_flag,T1,T2,context);
[Mx,My,Mz,dummy,dummy] = evolve_M_n_acq_2D_C(ParamsSt);
M1 = consolidate_3vec_Magnetization_to_2D_cell(Mx,My,Mz,Nx,Ny);
plot_2D_M_n_Phase(M0_sample,M1,'M(x,y) post RF','Phase post RF',context);

toc
return;

