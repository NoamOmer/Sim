
function [sample_z_axis,M_init] = UD1D_create_sample(context)
set_context;

declare_stage('Create sample');
% Create a new z-axis for the the sample
sample_z_axis = (sample_Zi : dz : sample_Zf+dz);

% % Gaussians
% loc = 400;
% % sample = exp(-70E+21*((sample_z_axis *dz).^8));
% sample = exp(-1.0E+1*((sample_z_axis+loc*dz - 50*dz).^2)) + ...
% 	     exp(-1.0E+1*((sample_z_axis-loc*dz + 250*dz).^2));% + ...
% % 		 exp(-2E+1*((sample_z_axis- 50*dz - 50*dz).^2));      % Smooth shape (gaussian) decreases the effect of Gibbs phenomenon

% Square
sample = exp(-1E+45*(((sample_z_axis - 1.4) * dz).^18)) + ...
         exp(-1E+45*(((sample_z_axis + 1.4) * dz).^18));
sample = sample .* exp(-1E+44*(((sample_z_axis) * dz).^18));

% Series of Squares
sample = 0;
NN = 8;
LL = 1.0*exp(-0.05*(1:NN)) .* linspace(sample_z_axis(1),sample_z_axis(end),NN);
for idx = 1:NN
	sample = sample + exp(-10^(61)*(((sample_z_axis + LL(idx)) * dz).^18));
end;
sample = sample .* exp(-1E+44*(((sample_z_axis) * dz).^18));
sample = transpose(smooth(sample,40));
% figure; plot(sample_z_axis,sample,'.-');

% % Triangle
% F = 500;
% sample = [zeros(1,F), 1:(length(z_axis)/2-F), (length(z_axis)/2-F):-1:1, zeros(1,F)]; sample = sample/max(sample);
% if (length(sample) < length(z_axis))   sample(end+1) = sample(end);   end;

% % Ones
% % sample = ones(1,length(sample_z_axis));

sample = 0.1 * sample ./ max(sample);
sample = sample + 1E-10;
M_init = transpose(sample) * M0;
if (DEBUG_FLAG >= 0)
	figure; plot(sample_z_axis,M_init(:,3),'.-'); title('Sample(z)'); xlabel('z-axis [cm]'); ylabel('M'); set_gca;
end;

return;
