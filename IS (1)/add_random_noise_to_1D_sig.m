% Input
% sig                    original 1D signal
% BW                     acquisition bandwidth (Hz)
% req_noise_std          required noise std in % of maximal signal value
% input_signal_lvl       reference signal level. Default - maximal value of signal
% complex_f              enable/disable complex noise
% reset_random_seed_f    enable/disable reset of random seed

function [sig_w_noise,noise_vec] = add_random_noise_to_1D_sig(sig,BW,req_noise_std,input_signal_lvl,complex_f,reset_random_seed_f,plot_f)

stream = RandStream.getDefaultStream;
if (reset_random_seed_f)
	stream.reset;
end;

if (input_signal_lvl > 0)
	sig_lvl = input_signal_lvl;
else
	sig_lvl = max(abs(sig));
end;

if (complex_f == 1)
	noise_vec = wgn(1,length(sig),sqrt(BW),1,stream,'complex');
else
	noise_vec = wgn(1,length(sig),sqrt(BW),1,stream,'real'   );
end;

if (1)

% Scale to required std level and signal value
noise_vec_std = sqrt((noise_vec * ctranspose(noise_vec)) / length(noise_vec));  % < n | n' >
noise_vec     = noise_vec * ((req_noise_std/100) * sig_lvl) / noise_vec_std;

else
	
noise_real = rand(1,length(sig));    % [ 0.0 ..  1.0]
noise_imag = rand(1,length(sig));

noise_real = noise_real - 0.5;       % [-0.5 .. +0.5]
noise_imag = noise_imag - 0.5;

noise_real_std = sqrt((noise_real * ctranspose(noise_real)) / length(noise_real));  % < n | n' >
noise_imag_std = sqrt((noise_imag * ctranspose(noise_imag)) / length(noise_imag));

% Remove mean value
noise_real = noise_real - mean(noise_real);
noise_imag = noise_imag - mean(noise_imag);

% Scale to required std level and signal value
noise_real = noise_real * ((req_noise_std/100) * sig_lvl) / noise_real_std;
noise_imag = noise_imag * ((req_noise_std/100) * sig_lvl) / noise_imag_std;

noise_vec = noise_real + 1i*noise_imag*complex_f;

end;
sig_w_noise = sig + noise_vec;

if plot_f
	figure; subplot(311); hold on;
	plot(abs(sig)          ,'r-');
% 	plot(abs(noise_vec)    ,'b-');
	plot(abs(sig+noise_vec),'k.');
	legend({'Sig','Sig+Noise'}); title('abs');
	subplot(312); hold on;
	plot(real(sig)         ,'k.-');
	plot(real(noise_vec)   ,'b.-');
	legend({'sig','noise'}); title('real');
	subplot(313); hold on;
	plot(imag(sig)         ,'k.-');
	plot(imag(noise_vec)   ,'b.-');
	legend({'sig','noise'}); title('imag');
end;

return

