% z_axis        : [cm]   : Excitation range
% Tp            : [sec]  : Theoretical excitation pulse duration
% rfwdth        : [sec]  : Actual excitation pulse duration
% Phi_e_sq      : [rad]  : Square part of post excitation phase
% OmegaE        : [Hz]   : omega(z) - required frequency response
% OmegaE_axis   : [Hz]   : Cartesian frequency axis
% ze_of_OmegaE  : [cm]   : z(omega) - inverted frequency responce
% swHz          : [Hz]   : Spectral width of the SLR pulse
% isPI          : [none] : Flag: is PI pulse (default is pi/2)
% slr_esp       : [none] : Edge slope parameter for the SLR pulse
 
function [pulse,rfpwr,rfpwrf] = convert_pulse_to_SLR(z_axis,Tp,rfwdth,Phi_e_sq,OmegaE,OmegaE_axis,...
                                                     ze_of_OmegaE,swHz,isPI,tpwr90,pw90,slr_esp,plot_flag)
set_globals;
swkHz = swHz*1E-3;

% dbstop in convert_pulse_to_SLR at 17
% Transform Phi(z)    --> Phi(z(w))
if (sum(ze_of_OmegaE - z_axis) < 1E-8)
    Phi_OmegaE = Phi_e_sq;
else
    Phi_OmegaE = interp1(z_axis,Phi_e_sq,ze_of_OmegaE);    % Only square part of the phase
end;

% --------------------------------------------------------------------------------------------
% The designed frequency domain filter (i.e. the pulse frequency response) is "casual". 
% This means that in the TIME domain the filter cannot begin before time zero.
% Therefore the time response starts at time zero, hence is not centered around zero.
% This shift in the time domain is analogue to a linear phase term in the frequency domain.
% In order to still get a pulse, which is symmetrical, we add a relevant linear shift in the
% required frequency responce:
% THIS IS DONE INSIDE THE SLR FUNCTION (INSIDE find_polyB)
% --------------------------------------------------------------------------------------------

% % Add re-focusing to the pulse
% if (RH_flag)
%     Phi_OmegaE = Phi_OmegaE + 2*pi*OmegaE_axis*Tp/2;
% else
%     Phi_OmegaE = Phi_OmegaE - 2*pi*OmegaE_axis*Tp/2;
% end;

% Spread Phi(z(w)) on a cartisian axis
w1 = linspace(OmegaE(1),OmegaE(end),length(OmegaE));
if (sum(OmegaE - w1) < 1E-8)
    Phi_omega = Phi_OmegaE;
else
    Phi_omega = interp1(OmegaE,Phi_OmegaE,w1);
end;

SRkHz     = abs(w1(end)-w1(1))*1E-3;
freq_resp = window_freq_resp(swkHz,SRkHz,Tp*1E+3,slr_esp,DEBUG_FLAG);
freq_resp = abs(freq_resp) + 1E-10;

if (1)
	% % % Create a smooth windowed frequency response [0..1]. Remove the phase,
	% % % imparted by the windowing function -- we need only the absolute value
	mid = round(length(freq_resp) /2);
	margin = 15;
	v     = freq_resp(mid-round(mid/margin) : mid + round(mid/margin));
	v_len = length(v);
	v     = v - 0.5*transpose(gausswin(v_len,4));
	freq_resp(mid-round(mid/margin) : mid + round(mid/margin)) = v;
	figure;plot(freq_resp,'.-');
end;

% Distribute Phi(z(w)) on a larger frequency range w/ size, equal to frequency response
dw = w1(2) - w1(1);
w2 = linspace(-sign(dw)*swHz/2, sign(dw)*(swHz/2 - abs(dw)), length(freq_resp));
Phi_omega = interp1(w1,Phi_omega,w2,'linear','extrap');
if (DEBUG_FLAG >= 3)
figure; plot(w1,Phi_OmegaE,'.-',w2,Phi_omega,'r--'); legend({'Phi_OmegaE(w1)','Phi_omega(w2)'});
end;

% Set the desired phase, on the frequency response amplitude 
freq_resp = freq_resp .* exp(i*Phi_omega);

freqkHz  = w2*1E-3;
if (DEBUG_FLAG >= 3)
	figure;
	subplot(3,1,1); plot(freqkHz,abs(freq_resp));
	title('Frequency Response (Amplitude)');  xlabel('Freq [kHz]'); ylabel('Amplitude'); set_gca;
	subplot(3,1,2); plot(freqkHz,phase(freq_resp));
	title('Frequency Response (Phase)');      xlabel('Freq [kHz]'); ylabel('Phase');     set_gca;
	subplot(3,1,3); plot(freqkHz,angle(freq_resp));
	title('Frequency Response (Angle)');      xlabel('Freq [kHz]'); ylabel('Phase');     set_gca;
end;

% For direct SLR interface, scale to freqeucncy axis to [-1..1]
freqkHz=2*freqkHz/swkHz;

% Note: pulse duration can (and should!) be longet than Tp
rfwdth_ms = rfwdth*1e+3;

dt_ms = 1/swkHz;
if (DEBUG_FLAG >= 3)
    figure;
    subplot(4,1,1); plot(freqkHz); title('freqkHz');
    subplot(4,1,2); plot(abs(freq_resp)); title('freq_resp');
    subplot(4,1,3); plot(phase(freq_resp)); title('phase freq_resp');
    subplot(4,1,4); plot(angle(freq_resp)); title('angle freq_resp');
end;

% Compute the SLR pulse (in units of [rad/ms]. Conversion to kHz can be done using divsion by 2*pi)
[pulse,A,B] = SLR(freqkHz,freq_resp,ones(1,length(freqkHz)),dt_ms,rfwdth_ms,isPI); % RF [rad/ms] ; dwell [ms] ; Tslr [ms]

% Since the SLR function is compatible with the 300MHz (which reverses the RF phase
% that is passed to it) we need to reverese the phase back to it's required value,
% when not passing the pulse to the 300MHz system.
if (SIMULATION)
    pulse = conj(pulse);
end;

z = exp(i*pi*freqkHz);

result1=polyval(flipdim(A,2),z.^(-1)); %The flipdim is due to the polyval function that takes the first term to multiply the highest power.
result2=polyval(flipdim(B,2),z.^(-1));
if (DEBUG_FLAG >= 3)

    if (isPI == 0)
        figure;
        subplot(2,1,1); plot(freqkHz*swkHz/2,2*abs(conj(result1).*result2));
        title('SLR: Frequency Response (Amplitude)');  xlabel('Freq [kHz]'); ylabel('Amplitude'); set_gca;

        phase1=phase(conj(result1).*result2)+(length(B)-1)/2*pi*freqkHz+phase(freq_resp);
        subplot(2,1,2); plot(freqkHz*swkHz/2,phase1-phase1(round(length(phase1)/2)));
        title('SLR: Frequency Response DIFF (Phase)'); xlabel('Freq [kHz]'); ylabel('Phase');     set_gca;
    else
        figure;
        subplot(2,1,1); plot(freqkHz*swkHz/2,abs(result2.^2));
        title('SLR: Frequency Response (Amplitude)');  xlabel('Freq [kHz]'); ylabel('Amplitude'); set_gca;

        phase1=phase(result2.^2)+2*(length(B)-1)/2*pi*freqkHz+phase(freq_resp);
        subplot(2,1,2); plot(freqkHz*swkHz/2,phase1-phase1(round(length(phase1)/2)));
        title('SLR: Frequency Response DIFF (Phase)'); xlabel('Freq [kHz]'); ylabel('Phase');     set_gca;
    end;
	figure;
	te_ms = linspace(0,rfwdth,length(pulse)) * 1E+3;
    subplot(3,1,1); plot(te_ms,  abs(pulse));  title('Pulse amplitude');  xlabel('Time [ms]');  set_gca;
	subplot(3,1,2); plot(te_ms,phase(pulse));  title('Pulse phase'    );  xlabel('Time [ms]');  set_gca;
	subplot(3,1,3); plot(te_ms,angle(pulse));  title('Pulse angle'    );  xlabel('Time [ms]');  set_gca;
end;

% Ugly
global FIDnum;
xd = freqkHz*swkHz/2;
yd = abs(result2.^2);
fn = sprintf('exc_amp_comp_%3.0f.mat',FIDnum);
save(fn,'xd','yd');

% Calibrate the pulse
cur_pwd = pwd;
cd 'D:\PhD\Matlab\PulseCalibration';
ampkHz = max(abs(pulse))/(2*pi);
[rfpwr,rfpwrf] = calib_power2(tpwr90,pw90,ampkHz);
rfpwrf = round(rfpwrf);
cd(cur_pwd);

% Convert the pulse to [T]
pulse = 1E+3 * pulse /gamma_T;  % [T] 1e+3: [rad/ms] --> [rad/sec]   ;   [rad/sec] / [rad/(sec*Tesla)] = [T]

return;

%----------
% % % Junk
%----------
% % % Reduce sample amplitude around phase non-continuity points
% mid = round(length(freq_resp) /2);
% freq_resp(mid-round(mid/200) : mid + round(mid/200)) = min(freq_resp);
% for smooth_idx = 1:3
% 	freq_resp = transpose(smooth(freq_resp,50));
% end;

