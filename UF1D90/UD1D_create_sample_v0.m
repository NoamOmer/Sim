
function [sample_z_axis,M_init] = UD1D_create_sample(context,SPACE_flag)
set_context;

declare_stage('Create sample');

% sample1 = zeros(1,length(sample_z_axis));

% % -----------------------------------------------------
% % % Gaussians
% % -----------------------------------------------------
% sample1 = exp(-0.1e-8*((sample_z_axis-1).^14));
% % sample1 = exp(-0.1e+1*((sample_z_axis-0).^16));
% 
% % sample1 = sample1 + 1 - 0.8*exp(-1e2*((sample_z_axis).^4));
% % figure; plot(sample1,'.-');
% loc = 10000;
% sample1 = exp(-1.0E+1*((sample_z_axis+loc*sample_dz).^2)) + ...
% 	      exp(-1.0E+1*((sample_z_axis-loc*sample_dz).^2));
% % figure; plot(sample1,'.-');
% % error('1');

% -----------------------------------------------------
% A SET OF GAUSSIANS WITH DECREASING WIDTHS
% -----------------------------------------------------
% sample1  = zeros(1,sample_Nz); %[zeros(1,0.2*Nz), ones(1,0.2*Nz), zeros(1,Nz - 0.4*Nz)];
% sample1(0.1*sample_Nz:0.25*sample_Nz) = 1;

% L_sigma = -sample_Lz/3;
% sample1  = exp(-((sample_z_axis)/(0.75*L_sigma)).^70);   % Thesis Square
% sample1  = exp(-((sample_z_axis)/(1.00*L_sigma)).^4  );  % Thesis wide Gaussian
% sample1  = exp(-((sample_z_axis)/(0.15*L_sigma)).^2  );  % Thesis thin Gaussian

% sft     = 0.0;
% sample1  = 0.0     + exp(-((sample_z_axis + sft + 5.0*sample_Lz/10)/(1.00*L_sigma)).^20);
% sample1  = sample1 + exp(-((sample_z_axis + sft + 0.0*sample_Lz/10)/(1.00*L_sigma)).^4 );
% sample1  = sample1 + exp(-((sample_z_axis + sft - 4.0*sample_Lz/10)/(0.25*L_sigma)).^2 );


% sample1 =           exp(-1e+7  *(((sample_z_axis +   sample_z_axis(end)/2) *sample_dz).^2));
% sample1 = sample1 + exp(-3e+7  *(((sample_z_axis +   sample_z_axis(end)/8) *sample_dz).^2));
% sample1 = sample1 + exp(-6e+7  *(((sample_z_axis -   sample_z_axis(end)/8) *sample_dz).^2));
% sample1 = sample1 + exp(-10e+7 *(((sample_z_axis -   sample_z_axis(end)/4) *sample_dz).^2));
% sample1 = sample1 + exp(-60e+7 *(((sample_z_axis - 2*sample_z_axis(end)/4) *sample_dz).^2));
% sample1 = sample1 + exp(-100e+7*(((sample_z_axis - 3*sample_z_axis(end)/4) *sample_dz).^2));
% figure; plot(sample1,'.-');

% -----------------------------------------------------
% Smooth square / Gaussian
% -----------------------------------------------------
% sample1 = exp(-1e-14*((sample_z_axis).^14));                                 % SAVE FOR SPEN SNR PAPER ANSWERS
% sample2 = 1*exp(-50e-2*((sample_z_axis-0.25*sample_z_axis(end)).^2));    % SAVE FOR SPEN SNR PAPER ANSWERS
% sample1 = sample1 - sample2;
if (exist('SPACE_flag','var') && (SPACE_flag == 1))
sample1 = exp(-10*((sample_z_axis).^10));     % SAVE FOR SPEN SNR PAPER ANSWERS
else
sample1 = exp(-0.5*((sample_z_axis).^6));     % SAVE FOR SPEN SNR PAPER ANSWERS
end;
% figure; plot(sample1);

% % % % % SAVE FOR SPEN SNR PAPER ANSWERS - Compare @-Nyquist vs. non-Nyquist SR, SR wOnes; SE
% % % % sample1 = exp(-0.01*((sample_z_axis-0.00*sample_z_axis(end)).^4));
% % % % sample2 = exp(-0.4*((sample_z_axis-0.55*sample_z_axis(end)).^4));
% % % % sample3 = exp(-0.4*((sample_z_axis+0.65*sample_z_axis(end)).^4));
% % % % sample1 = sample1 + 0.5*sample2 + 0.25*sample3;

% % -----------------------------------------------------
% % Sharp square
% % -----------------------------------------------------
% sample1 = zeros(1,length(sample_z_axis));
% mid  = round(length(sample_z_axis)/2);
% sample1(mid) = 1;
% span = 5000;
% sample1([mid-span,mid+span]) = 1;
% sample1([mid-span:mid+span]) = 1;
% % figure; plot(sample1,'.-');

% -----------------------------------------------------
% Point
% -----------------------------------------------------
% Two points
% loc = [2000,1500,1000,500,400,350,300];
% loc = loc(5);
% sample1(round(length(sample_z_axis)/2)+loc) = 1;
% sample1(round(length(sample_z_axis)/2)-loc) = 1;
% figure; plot(sample1,'s');

% % Single point
% sample1 = zeros(1,length(sample_z_axis));
% sample1(round(length(sample_z_axis)/2)) = 1;
% 
% -----------------------------------------------------
% Single pixel
% -----------------------------------------------------
% n_px_pts = length(sample_z_axis) * required_res / Lhalf*2;
% sample1 = zeros(1,length(sample_z_axis));
% mid     = round(length(sample_z_axis)/2);
% sample1(round(mid-n_px_pts/2) : round(mid+n_px_pts/2)) = 1;
% % figure; plot(sample1,'.-');

% -----------------------------------------------------
% Multiple squares
% -----------------------------------------------------
% sample1 = 0;
% NN = 8;
% LL = exp(0.1*(NN:-1:1)) .* linspace(sample_Zi+0.35*sample_Lz,sample_Zf-0.25*sample_Lz,NN);
% for idx = 1:NN
% 	sample1 = sample1 + exp(-10^(50+idx*16)*(((sample_z_axis + LL(idx)) * sample_dz).^(16+idx*2)));
% end;
% sample1 = transpose(smooth(sample1,40));
% figure; plot(sample1,'.');

% -----------------------------------------------------
% % Series of Squares
% -----------------------------------------------------
% sample1 = 0;
% NN = 6;
% LL = 1.0*exp(-0.05*(1:NN)) .* linspace(sample_z_axis(1),sample_z_axis(end),NN);
% for idx = 1:NN
% 	sample1 = sample1 + exp(-10^(61)*(((sample_z_axis + LL(idx)) * sample_dz).^14));
% end;
% sample1 = sample1 .* exp(-1E+44*(((sample_z_axis) * sample_dz).^14));
% sample1 = transpose(smooth(sample1,40));
% % figure; plot(sample_z_axis,sample1,'.-');

% -----------------------------------------------------
% % Triangle
% -----------------------------------------------------
% F = 500;
% sample1 = [zeros(1,F), 1:(length(sample_z_axis)/2-F), (length(sample_z_axis)/2-F):-1:1, zeros(1,F)]; sample1 = sample1/max(sample1);
% if (length(sample1) < length(sample_z_axis))   sample1(end+1) = sample1(end);   end;

% -----------------------------------------------------
% % Ones
% -----------------------------------------------------
% sample1 = ones(1,length(sample_z_axis));

% sample1 = 0.1 * sample1 ./ max(sample1);
sample1 = exp(-10*((sample_z_axis).^10));     % SAVE FOR SPEN SNR PAPER ANSWERS
sample1 = sample1 + 1E-6;
M_init = transpose(sample1) * M0;
% if (DEBUG_FLAG >= 2)
% if (exist('SPACE_flag','var') && (SPACE_flag == 1))
% 	figure; plot(sample_z_axis,M_init(:,3),'.-'); title('Sample(z)'); xlabel('z-axis [cm]'); ylabel('M'); set_gca;
% 	uiwait(msgbox('Please approve sample'));
% end;
% end;

return;
