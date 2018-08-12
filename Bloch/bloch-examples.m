% -----------------------------
%  G A R B A G E  --  bloch.m
% -----------------------------

Free Precession
-----------------
%duration = 10E-6;
%rad_freq = gamma_T*B0_lab(3);
%unit_vec = [0 0 1];
%M_rot = free_precess(M_rot, duration, dt, T1, T2, unit_vec, rad_freq);

Initial pulse application (using cross product)
-----------------------------------------------
% 	for idx = 1:length(t)
% 		dM(idx,:)  = dt * gamma_T * cross(M_rot(idx,:), B_rot(idx,:));
% 		M_rot(idx+1,:) = M_rot(idx,:) + dM(idx,:);
% 
% 		M_rot(idx+1,:) = T1_relaxation(M_rot(idx+1,:), 1/2, dt, T1);
% 		M_rot(idx+1,:) = T2_relaxation(M_rot(idx+1,:), dt, T2);
% 	end;
% 	plot([0 t],M_rot(:,1),'r.-', [0 t],M_rot(:,2),'b.-', [0 t],M_rot(:,3),'k.-');
% 	legend({'Mx\_rot','My\_rot','Mz\_rot'});


Initial application of pulse using explicit form of bloch equations
-------------------------------------------------------------------
% M0 = (1/2)*[0,0,1];            % [amper/m] Magnetization
% Mx_rot(1) = M0(1);
% My_rot(1) = M0(2);
% Mz_rot(1) = M0(3);
% for idx = 1:length(t)
% 	dMx(idx) = dt * gamma_T * ( My_rot(idx)*Bz_rot(idx) - Mz_rot(idx)*By_rot(idx) );
% 	dMy(idx) = dt * gamma_T * ( Mz_rot(idx)*Bx_rot(idx) - Mx_rot(idx)*Bz_rot(idx) );
% 	dMz(idx) = dt * gamma_T * ( Mx_rot(idx)*By_rot(idx) - My_rot(idx)*Bx_rot(idx) );
% 
% 	Mx_rot(idx+1) = Mx_rot(idx) + dMx(idx);
% 	My_rot(idx+1) = My_rot(idx) + dMy(idx);
% 	Mz_rot(idx+1) = Mz_rot(idx) + dMz(idx);
% end;
% figure; plot([0 t],Mx_rot,'r-.');
% figure; plot([0 t],My_rot,'b--');
% figure; plot([0 t],Mz_rot,'k.-');

% plot([0 t],Mx_rot,'r-.',[0 t],My_rot,'b--',[0 t],Mz_rot,'k.-');
% legend({'Mx_rot','My_rot','Mz_rot'});

