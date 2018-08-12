% Input
% sig                    original 1D signal
% BW                     acquisition bandwidth (Hz)
% req_noise_std          required noise std in % of maximal signal value
% input_signal_lvl       reference signal level. Default - maximal value of signal
% complex_f              enable/disable complex noise
% reset_random_seed_f    enable/disable reset of random seed

function [mat_w_noise,noise_mat] = add_random_noise_to_2D_mat(mat,BW,req_noise_std,input_signal_lvl,complex_f,reset_random_seed_f,plot_f)

stream = RandStream.getDefaultStream;
if (reset_random_seed_f)
	stream.reset;
end;

if (input_signal_lvl > 0)
	sig_lvl = input_signal_lvl;
else
	sig_lvl = max(abs(mat(:)));
end;

if (complex_f == 1)
	noise_mat = wgn(size(mat,1),size(mat,2),sqrt(BW),1,stream,'complex');
else
	noise_mat = wgn(size(mat,1),size(mat,2),sqrt(BW),1,stream,'real'   );
end;

% Scale to required std level and signal value
noise_mat_std = std(noise_mat(:));
noise_mat     = noise_mat * ((req_noise_std/100) * sig_lvl) / noise_mat_std;

mat_w_noise = mat + noise_mat;

if plot_f
	figure; subplot(311); hold on;
	plot(abs(mat)          ,'r-');
% 	plot(abs(noise_mat)    ,'b-');
	plot(abs(mat+noise_mat),'k.');
	legend({'Sig','Sig+Noise'}); title('abs');
	subplot(312); hold on;
	plot(real(mat)         ,'k.-');
	plot(real(noise_mat)   ,'b.-');
	legend({'sig','noise'}); title('real');
	subplot(313); hold on;
	plot(imag(mat)         ,'k.-');
	plot(imag(noise_mat)   ,'b.-');
	legend({'sig','noise'}); title('imag');
end;

return

