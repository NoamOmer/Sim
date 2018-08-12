
function [z_axis,P_DNu0,fh] = UD1D_retrieve_inhomogeneity(context)
set_context;
z_axis = 0;
fh     = 0;
inhomo_OP = 4;
switch inhomo_OP
	case 0
		P_DNu0 = 0;
	case 1
		P_DNu0 = [0 0  0 1]*2000; % constant off-resonance (value for paper [0 0 1]*2000)
	case 2
		P_DNu0 = [0 0  1 0]*100;  % Linear
	case 3
		P_DNu0 = [0 1  0 0]*015;  % ^2
	case 4
		P_DNu0 = [1 5 15 0]*00;  % ^3
end;

% declare_stage('Retrieve the B0 inhomogeneity and set the spatial axis');
% fh = 0;
% if (Map_Mult_Factor ~= 0)
% 	if (Seq1D)
% 		[map_fn,P_DNu0,z_axis,tof,op,mean_std] = Inhomo_Zmap;
% 	else
% 		map_fn = '31mar08_1_MATOrient_RP90_SF2_DSF1_PO5.mat';%uigetfile('*.mat','Select an Oriented Phase-Map File');
% 		load(map_fn);
% 		[z_axis,P_DNu0,Porder] = calc_oriented_1D_inhomo(map_mat, map_mask, x_axis, y_axis, rot_phi, 0, 0, gPolyOrder);
% 	end;
% 	if (DEBUG_FLAG >= 1)
% 		fh = figure; plot(z_axis,polyval(P_DNu0,z_axis)); hold on;
% 		title('Field inhomogeneity'); xlabel('z-axis [cm]'); ylabel('\Delta\Nu [Hz]'); set_gca;
% 	end;
% 	% Inhomogeneous field tests
% 	P_DNu0 = Map_Mult_Factor*P_DNu0;
% 	plot(z_axis,polyval(P_DNu0,z_axis),'c');
% else
% 	P_DNu0 = 0;
% 	z_axis = 0;
% end;

return;

