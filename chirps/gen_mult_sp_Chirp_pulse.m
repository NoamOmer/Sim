% Generate multiple stationary point chirp pulse.
% Uses the RH rotation convention.

% Input parameters
% te          : [sec]    : Excitation time vector (0 .. Tp)
% Oi          : [Hz]     : Initial excitation frequency
% R           : [Hz/sec] : Chirp rate - equal for all sub-chirps
%             :          : Should not be absolute-valued albeit contain a sign, signifying the sweeping direction
% Nsp         : [none]   : Number of stationary points
% B1ampCoeff  : [none]   : Coefficient for RF field strength
% Ge_         : [G/cm]   : Excitation gradient


function [chirp_rot,alpha0_n,alpha1_n,alpha2] = gen_mult_sp_Chirp_pulse(te,Oi,R,Nsp,B1ampCoeff,Ge_)
set_globals;

% if (Nsp == 1)
%     [chirp_rot,phi_chirp_rot,alpha0_n,alpha1_n,alpha2] = gen_Chirp_pulse(te,Oi,R,context);
%     return;
% end;

% Calculate the excited frequency range
Tp = te(end);                                       % [sec]
dO = R*Tp;                                          % [Hz]

% Add overlapp between pulses
overlap_factor = 0/100;                             % [%]
Oi     = Oi - overlap_factor * dO;                  % [Hz]
cur_R  = R  * (1 + 2*overlap_factor);               % [Hz/sec]
cur_dO = cur_R*Tp;                                  % [Hz]
B1ampLabHz  = B1ampCoeff*sqrt(abs(cur_R));          % [Hz]

% Calculate the pulses' amplitude
amp_chirp_rot = (1/2)*(2*pi*B1ampLabHz)/gamma_T;    % [T]

chirp_rot = 0;
% a = [];
for idx = 1:Nsp
	% Calculate the pulse frequency range
	Oi_n(idx) = Oi + (idx-1)*dO;                                         % [Hz]
	Of = Oi_n(idx) + cur_dO;                                             % [Hz]
    
	phi_chirp_rot = 2*pi * (Oi_n(idx)*te + (1/2)*cur_R*(te.^2));         % [rad]
    if (RH_flag)
    	chirp_rot = chirp_rot + amp_chirp_rot * exp(+1i*phi_chirp_rot);   % '+j' <--> Right-Hand rotated pulse
%       a(idx,:) = chirp_rot;
    else
    	chirp_rot = chirp_rot + amp_chirp_rot * exp(-1i*phi_chirp_rot);   % '-j' <--> Left -Hand rotated pulse
    end;
	% Calculate each sub-pulse post-excitation phase coefficients
    alpha0_n(idx) = - Tp * ((Oi_n(idx))^2) / (2*cur_dO);
	alpha1_n(idx) = gammaHz * Ge_ * Tp * Of / cur_dO;
    
    % Verify that in case of a single stationaty point we get the expected (reduced expression) values
    if (Nsp == 1)
        if ( (abs(alpha0_n(idx) - (- ((Oi_n(idx))^2) / (2*cur_R)))  > abs(alpha0_n(idx)/1000)) || ...
             (abs(alpha1_n(idx)  - (gammaHz * Ge_ * Tp * Of / cur_dO)) > abs(alpha1_n(idx)/1000)) )
            uiwait(errordlg('gen_mult_sp_Chirp_pulse: Inconsistent phase profile coefficients value'));
        end;
    end;
end;
% chirp_rot = a(1,:);

% The post-excitation quadratic phase coefficient is similar for all sub pulses,
% thus calculated outside the for-loop
alpha2 = - ((gammaHz*Ge_)^2) / (2*cur_R);

return;


