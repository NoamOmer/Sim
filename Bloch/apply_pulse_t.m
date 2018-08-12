
% Apply a magnetic field pulse on a single value magnetization vector.
% Final M is given through time.
function [M] = apply_pulse_t(B_eff_rot, M);

Bloch_set_globals;

%dM = 0;
for idx = 1:length(B_eff_rot)
	% 1st method - using simple form of Bloch equations.
	% Problematic method: causes a gradual increase of Mx & My while precessing in XY plane
	% (for long time durations T). It is a result of not using an infinitesimal
	% value of dt -- i.e., the smaller the dt --> the less significant is the problem
% 		dM(end,1:3)  = dt * gamma_T * cross(M_rot(end,:), B_rot(idx,:));
% 		M_rot(end+1,:) = M_rot(end,:) + dM(end,:);

	% 2nd method - using solution of Bloch equations, i.e., rotation around effective field.
    % (This method does not have the problem described for the 1st method)
	B_vec = B_eff_rot(idx,:);
	B_vec_sz = sqrt(sum(B_vec.^2));
    if (B_vec_sz == 0)
        R = [1 0 0; 0 1 0; 0 0 1];
    else;
		B_unit_vec = B_vec / B_vec_sz;
		phi_rot = gamma_T*B_vec_sz*dt;
		R = rot(B_unit_vec,-phi_rot);
	end;
	M(end+1,:) = (R*M(end,:)')';    %'
    
    % T1 & T2 relaxation
    if (Relax_Flag)
		M(end,:) = T1_relaxation(M(end,:), M0_rot(3), dt, T1);
		M(end,:) = T2_relaxation(M(end,:), dt, T2);
	end;
end;

return;

