% evolve_M_CPP.c
% evolve_M_n_acq1Dgrad_t_CPP.c
% evolve_M_n_acq_2D_C.c
% evolve_M_n_acq_2D_tIndependent_C.c
% evolve_M_n_acq_wTimeDepGrad_2D_C.c

function compile_C_code()
set_globals;

cur_cd = pwd;

% ---------------------------------------------------------------------------------
% Compile I
% ---------------------------------------------------------------------------------
global SimRoot;
cd([SimRoot 'C code']);

filenames = {'evolve_M_CPP'                    ,...
             'evolve_M_n_acq1Dgrad_t_CPP'      };
%              'evolve_M_n_acq_2D_C'             ,...
%              'evolve_M_n_acq_2D_tIndependent_C',...
%              'evolve_M_n_acq_wTimeDepGrad_2D_C'};
header_fn = 'evolve.h';

for idx = 1:length(filenames)
	full_file_nm   = [SimRoot 'C code' filesep  filenames{idx} '.c'      ];
	full_dll_nm    = [SimRoot 'C code' filesep  filenames{idx} '.' mexext];
	full_header_fn = [SimRoot 'C code' filesep  header_fn                ];

	if (exist(full_dll_nm,'file') == 3)
		file_prop = dir(full_file_nm  );   file_date = datenum(file_prop.date);
		dll_prop  = dir(full_dll_nm );     dll_date  = datenum(dll_prop.date);
		hdr_prop  = dir(full_header_fn);   hdr_date  = datenum(hdr_prop.date);
		if (file_date < dll_date) && (hdr_date < dll_date)
			continue;
		end;
	end;
	
	fprintf('compiling %s\n',full_file_nm);
	mex(full_file_nm);
	
	
% 	mex('-llibemlrt'                ,...
% 	    '-llibeng'                  ,...
% 	    '-llibfixedpoint'           ,...
% 	    '-llibmat'                  ,...
% 	    '-llibmex'                  ,...
% 	    '-llibmwblas'               ,...
% 	    '-llibmwblascompat32'       ,...
% 	    '-llibmwcgir_construct'     ,...
% 	    '-llibmwlapack'             ,...
% 	    '-llibmwmathutil'           ,...
% 	    '-llibmwservices'           ,...
% 	    '-llibmwsl_fileio'          ,...
% 	    '-llibmwsl_solver_rtw'      ,...
% 	    '-llibmx'                   ,...
% 	    '-llibut'                   ,...
% 	    '-lmclbase'                 ,...
% 	    '-lmclmcr'                  ,...
% 	    '-lmclmcrrt'                ,...
% 	    '-lne_mli'                  ,...
% 	    '-lne_rtl'                  ,...
% 	    '-lphysmod_foundation_util' ,...
% 	    '-lrtwcg'                   ,...
% 		'-lliba'                    ,...
% 	    '-LC:\Program Files\MATLAB\R2010b\extern\lib\win64\microsoft',...
% 	    '-LE:\01 Post\Matlab\ChirpGenerator\cSLR2'                   ,...
% 		full_file_nm);
	
end;

% ---------------------------------------------------------------------------------
% Compile II
% ---------------------------------------------------------------------------------
global EMC_root;
fn = {[EMC_root 'Algorithm' filesep 'C_code' filesep 'emc_fit_kernel_C']};
for idx = 1:length(fn)
	filename = fn{idx};
	full_file_nm  = [filename '.c'          ];

	full_dll_nm = [filename '.' mexext];	
	
	compile_flag = 0;
	if (exist(full_dll_nm,'file'))
		file_prop = dir(full_file_nm);    file_date = datenum(file_prop.date);
		dll_prop  = dir(full_dll_nm);     dll_date  = datenum(dll_prop.date);
		if (file_date >= dll_date)
			compile_flag = 1;
		end;
	else
		compile_flag = 1;
	end;
	if (compile_flag)
		fprintf('compiling %s\n',full_file_nm);
		mex(full_file_nm);
	end;
end;

% ---------------------------------------------------------------------------------
% done
% ---------------------------------------------------------------------------------
cd(cur_cd);

return;


