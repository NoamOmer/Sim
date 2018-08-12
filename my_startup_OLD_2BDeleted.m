disp('------------------------------------------------------------------')
disp('                        WELCOME TO MATLAB                         ')
disp('------------------------------------------------------------------')
disp(' ')
disp('Practice does not make perfect. Only perfect practice makes perfect.')
disp(' ')
disp(' ')
disp(' ')
disp('what''s on your mind?')
disp(' ')
disp(' ')
disp(' ')
% pause(0);

% ===================
%  DEBUG OPTIONS
% ===================
dbclear all;

% ===================
%  ADD Path
% ===================
cd('/Users/noambe/Dropbox/1/Sim/');
set_globals;
cd('/Users/noambe/Dropbox');
Hybrid_Siemens = 1;


addpath([Postroot 'Matlab' filesep 'ChirpGenerator'                 ]);
addpath([Postroot 'Matlab' filesep 'ChirpGenerator' filesep 'cSLR' ]);
addpath([Postroot 'Matlab' filesep 'ChirpGenerator' filesep 'cSLR2']);

addpath([Postroot 'Matlab' filesep 'SPM' filesep 'spm8']);
% addpath([Postroot 'Matlab' filesep 'SPM' filesep 'NIFTI_20110921']);

addpath([SimRoot                          ]);
addpath([SimRoot '2D pulses'              ]);
addpath([SimRoot '2dHybridUFEPI'          ]);
addpath([SimRoot 'Bloch'                  ]);
addpath([SimRoot 'C code'                 ]);
addpath([SimRoot 'chirps'                 ]);
addpath([SimRoot 'COMMON'                 ]);
addpath([SimRoot 'DICOMviewer'            ]);
addpath([SimRoot 'EPI'                    ]);
addpath([SimRoot 'GE'                     ]);
addpath([SimRoot 'IS'                     ]);
addpath([SimRoot 'IS' filesep 'backcor'   ]);
addpath([SimRoot 'IS' filesep 'NIFTI_API' ]);
addpath([SimRoot 'Magnetization Evolution']);
addpath([SimRoot 'MultSP'                 ]);
addpath([SimRoot 'RF Pulse'               ]);
addpath([SimRoot 'RAISE'                  ]);
addpath([SimRoot 'SPA'                    ]);
addpath([SimRoot 'SR'                     ]);
addpath([SimRoot 'T2compHybridUFEPI'      ]);
addpath([SimRoot 'TSE'                    ]);
addpath([SimRoot 'UF1D90'                 ]);
addpath([SimRoot 'UF1D180'                ]);
addpath([SimRoot 'Siemens'                ]);
addpath([SimRoot 'DICOMviewer'            ]);
addpath([SimRoot 'OMP'                    ]);

addpath([ReadSiemensDir                   ]);


addpath([Projectsroot 'Pulse_Design']);
addpath([Projectsroot 'Siemens' filesep 'read .dat files']);
addpath([Projectsroot 'Siemens' filesep 'read .dat files' filesep 'Jons_Functions']);
addpath([Projectsroot '2DEPI']);
addpath([Projectsroot 'Diffusion']);
addpath([Projectsroot 'CS EPI']);
addpath([Projectsroot 'CS EPI' filesep 'cs_sampling_mask']);
addpath([Projectsroot 'CS EPI' filesep 'Undersampling masks']);
addpath([Projectsroot 'CS EPI' filesep 'demo_cs_mri']);
addpath([Projectsroot 'CS EPI' filesep 'spatial_fft']);
addpath([Projectsroot 'CS EPI' filesep 'EPI_reco']);
addpath([Projectsroot 'CS EPI' filesep 'sparseMRI_v0.2']);
addpath([Projectsroot 'CS EPI' filesep 'sparseMRI_v0.2' filesep 'utils']);
addpath([Projectsroot 'SPEN SNR']);
addpath([Projectsroot 'Prostate']);
addpath([Projectsroot 'Prostate' filesep 'CS_Recon']);
addpath([Projectsroot 'Prostate' filesep 'CS_Recon' filesep '@MCNUFFT']);
addpath([Projectsroot 'Prostate' filesep 'CS_Recon' filesep '@TempPCA']);
addpath([Projectsroot 'Prostate' filesep 'CS_Recon' filesep 'nufft_files']);
addpath([Projectsroot 'Prostate' filesep 'CS_Recon' filesep 'utils']);

addpath('/Users/noambe/Dropbox/EMC_T2_FIT/Algorithm/');
addpath('/Users/noambe/Dropbox/1/MS_project/');
addpath('/Applications/freesurfer/matlab/');
addpath('/Users/noambe/Dropbox/1/Multi_Comp_T2/');
addpath('/Users/noambe/Dropbox/1/MultiCompT2_AmirBeck/');

if (~isMac)
	% Wavelets
	run('C:\Program Files\MATLAB\R2010b\toolbox\Wavelab850\WavePath.m');
end;

if (Hybrid_Siemens)
	addpath([Projectsroot 'Siemens' filesep]);

	addpath([SiemensRoot              ]);
	addpath([SiemensRoot 'GE_reco'    ]);
	addpath([SiemensRoot 'EPI_reco'   ]);

	rmpath([Projectsroot '2DHybridUF' ]);
else
	rmpath([SiemensRoot               ]);
	rmpath([SiemensRoot 'GE_reco'     ]);
	rmpath([SiemensRoot 'EPI_reco'    ]);

	addpath([Projectsroot '2DHybridUF']);
end;


return;
% ===================
%  Screen properties
% ===================
% - - - - - - - - -
% In native units:
% - - - - - - - - -
% set(0,'defaultfigureposition', [2   65   1025   630]);

% - - - - - - - - - - - - -
% Or, in normalized units:
% - - - - - - - - - - - - -
% set(0,'defaultfigureunits'   , 'norm');
% set(0,'defaultfigureposition', [.002 0.058 .997 .87]);

%%
set(0,'defaultfigureunits','points');   % for MATLAB packages that opens figures in strange screen locations
% set(0,'defaultfigureposition', [200 200 1035 600]);
f=200;
set(0,'defaultfigureposition', [1+f 1+f 1385-f 950-f]); % full frame on mac
% figure;
%%

% - - - - - - - - - - - - - - - - - - - - - -
% Printing preferences - for copying to word
% - - - - - - - - - - - - - - - - - - - - - -
set(0,'DefaultFigurePaperType','A4');
set(0,'DefaultFigurePaperOrientation','portrait');
set(0,'DefaultFigurePaperUnits','inches');
set(0,'DefaultFigurePaperPosition',[0 0 5 3]);

% % Printing preferences - for A4 Landscape, full scale
% % - - - - - - - - - - - - - - - - - - - - - - - - - -
% set(0,'DefaultFigurePaperType','A4');
% set(0,'DefaultFigurePaperOrientation','landscape');
% set(0,'DefaultFigurePaperUnits','inches');
% set(0,'DefaultFigurePaperPosition',[0.2500 0.4612 11.1929 7.3453]);

