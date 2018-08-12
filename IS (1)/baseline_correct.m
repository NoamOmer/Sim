function [fid_mat_both_bl_corrected] = baseline_correct(fid_mat_both,bl_ord)

bl_margin_precentage = 0;

v_abs =  abs(fid_mat_both(round(size(fid_mat_both,1)/2),:));
v_rl  = real(fid_mat_both(round(size(fid_mat_both,1)/2),:));
v_im  = imag(fid_mat_both(round(size(fid_mat_both,1)/2),:)); % figure; subplot(311); plot(v_abs); subplot(312); plot(v_rl); subplot(313); plot(v_im);
[max_val_abs,max_loc_abs] = max(v_abs);
[max_val_rl ,max_loc_rl ] = max(v_rl );
[max_val_im ,max_loc_im ] = max(v_im );

bl_margin    = round( size(fid_mat_both,2) * bl_margin_precentage / 100);
bl_idx_abs   = [1 : (max_loc_abs-bl_margin), (max_loc_abs+bl_margin) : size(fid_mat_both,2)];
bl_idx_rl    = [1 : (max_loc_rl -bl_margin), (max_loc_rl +bl_margin) : size(fid_mat_both,2)];
bl_idx_im    = [1 : (max_loc_im -bl_margin), (max_loc_im +bl_margin) : size(fid_mat_both,2)];
bl_trshd_abs = max(v_abs(bl_idx_abs));
bl_trshd_rl  = max(v_rl (bl_idx_rl ));
bl_trshd_im  = max(v_im (bl_idx_im ));
N            = length(v_abs);
% bl_n         = linspace(-round(N/2),+round(N/2),N);
bl_n         = linspace(1,+round(N),N);
% bl_ord       = 0;
bl_fct       = 'stq';   % [sh | ah | stq | atq]

for idx = 1:size(fid_mat_both,1)
	v1     =      fid_mat_both(idx,:);
	v1_abs =  abs(fid_mat_both(idx,:));
	v1_rl  = real(fid_mat_both(idx,:));  % v1_rl_mean = mean(v1_rl);    v1_rl = v1_rl - v1_rl_mean;
	v1_im  = imag(fid_mat_both(idx,:));  % v1_im_mean = mean(v1_im);    v1_im = v1_im - v1_im_mean;
	
	if (idx == 22)
		% 		figure; subplot(211); plot(abs(fftshift(fft(v1)))); subplot(212); plot(abs(fftshift(fft(v1 - mean(v1)))),'r');
	end;
	% complex
	[bl,a,it,ord,s,fct] = backcor(bl_n,v1,bl_ord,bl_trshd_abs,bl_fct);
	bl     = transpose(bl);
	% 	bl     = bl - mean(bl);
	v2     = v1 - bl;            % figure; subplot(211); plot(v1_abs,'b'); hold on; plot(bl_abs,'r'); subplot(212); plot(abs(v2_abs),'b');
	
	% abs
	[bl_abs,a,it,ord,s,fct] = backcor(bl_n,v1_abs,bl_ord,bl_trshd_abs,bl_fct);
	bl_abs = transpose(bl_abs);
	% 	bl_abs = bl_abs - mean(bl_abs);
	v2_abs = v1_abs - bl_abs;
	v2_abs = v2_abs .* exp(1i*angle(v1));   % figure; subplot(211); plot(v1_abs,'b'); hold on; plot(bl_abs,'r'); subplot(212); plot(abs(v2_abs),'b');
	
	% real
	[bl_rl,a,it,ord,s,fct] = backcor(bl_n,v1_rl,bl_ord,bl_trshd_rl,bl_fct);
	bl_rl = transpose(bl_rl);
	% 	bl_rl = bl_rl - mean(bl_rl);
	v2_rl = v1_rl - bl_rl/1;  % figure; subplot(211); plot(v1_rl,'b'); hold on; plot(bl_rl,'r'); subplot(212); plot(v2_rl,'b');
	
	% imag
	[bl_im,a,it,ord,s,fct] = backcor(bl_n,v1_im,bl_ord,bl_trshd_im,bl_fct);
	bl_im = transpose(bl_im);
	% 	bl_im = bl_im - mean(bl_im);
	v2_im = v1_im - bl_im/1;   %figure; subplot(211); plot(v1_im,'b'); hold on; plot(bl_im,'r'); subplot(212); plot(v2_im,'b');
	
% 	fid_mat_both_bl_corrected(idx,:) = v1 - mean(v1);
		fid_mat_both_bl_corrected(idx,:) = v2;
	% 	fid_mat_both_bl_corrected(idx,:) = v2_abs;
	% 	fid_mat_both_bl_corrected(idx,:) = (v2_rl + v1_rl_mean) + 1i*(v2_im + v1_im_mean);
	% 	fid_mat_both_bl_corrected(idx,:) = v2_rl + 1i*v2_im;
end;

return;

