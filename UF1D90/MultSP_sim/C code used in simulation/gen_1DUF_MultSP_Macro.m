
function gen_1DUF_MultSP_Macro(void)
clcl;
set_gen_1DUF_MultSP_Macro_globals;

% gammaHz = 4.2574e+3;
% sw = 250e+3;
% Tp = 3e-3;
% Ta = Tp;
% Ge = 7.8;
% dta = 1/sw;
% Lz  = 7.5;

nMz = Ta/dta;
Dk = 2*nMz/Lz;
k0_tag = gammaHz*Ge*Tp;
disp(sprintf(' Dk = %3.3f\n k0_tag = %3.3f\n',Dk,k0_tag));
% return;
% --------------------------------------------------------------------------
% Create file names
% --------------------------------------------------------------------------
cc = fix(clock);
dd = num2str(cc(3));
mm = num2str(cc(2));
yy = num2str(cc(1));  yy = yy(3:4);
cc = [dd,mm,yy,'_',num2str(cc(4)),num2str(cc(5)),num2str(cc(6))];

chirp_name  = sprintf('nbe_chrp90_%s.RF',cc);
macro_name  = sprintf('nbe_macro_%s'    ,cc);

chirp_fn = [root_dir chirp_name];
macro_fn = [root_dir macro_name];

% --------------------------------------------------------------------------
% RF
% --------------------------------------------------------------------------
declare_stage('Design RF pulse');
OmegaE      = gammaHz*Ge*z_axis;
OmegaE_axis = OmegaE;
z_of_OmegaE = interp1(OmegaE,z_axis,OmegaE);
swHz        = 1/dte;
isPI        = 0;
alpha2 = - ((gammaHz*Ge)^2) / (2*rf_R);
if (Nsp > 1)
	n = round(linspace(0,length(z_axis),Nsp+1));
	Phi_e_sq_rad = [];
	for sp_idx = 1:Nsp
		ax = z_axis(n(sp_idx)+1 : n(sp_idx+1));                               % Spatial range of current SP
		mid_val = mean([z_axis(n(sp_idx)+1), z_axis(n(sp_idx+1))]);           % Middle value of the range
		Phi_e_sq_rad(end+1:n(sp_idx+1)) = 2*pi*alpha2*((ax - mid_val) .^ 2);  % Frequency response
	end;
else
	Phi_e_sq_rad = 2*pi*alpha2*(z_axis.^2);
end;

if (DEBUG_FLAG)
figure;
subplot(3,1,1); plot(z_axis     ,OmegaE*1E-3 ,'.-'); title('OmegaE(z)  '); xlabel('z-axis [cm]'   ); ylabel('Omega(z) [kHz]'       ); grid;
subplot(3,1,2); plot(OmegaE*1E-3,z_of_OmegaE ,'.-'); title('z(OmegaE)  '); xlabel('Omega(z) [kHz]'); ylabel('z-axis [cm]'          ); grid;
subplot(3,1,3); plot(z_axis     ,Phi_e_sq_rad,'.-'); title('Phi e sq(z)'); xlabel('z-axis [cm]'   ); ylabel('Phi_e_sq_rad(z) [rad]'); grid;
end;
declare_substage('SLR');
[chirp_rot,rfpwr,rffpwr] = convert_pulse_to_SLR(z_axis,Tp,rfwdth,Phi_e_sq_rad,OmegaE,OmegaE_axis,...
                                                z_of_OmegaE,swHz,isPI,tpwr90,pw90,RF_ESP_SLR,DEBUG_FLAG);

declare_substage('RF pulse file');
disp(sprintf('Saving RF pulse file (%s)',macro_fn));
create_chirp_table_varian(abs(chirp_rot),RFPhaseSign *phase(chirp_rot),chirp_fn);

% --------------------------------------------------------------------------
% Macro
% --------------------------------------------------------------------------
declare_stage('Macro file');
disp(sprintf('Generating MACRO file (%s)',macro_fn));
real_params_list = {'Nsp'            ,...
	                'rfwdth'         ,...
                    'rfpwr'          ,...
                    'rffpwr'         ,...
                    'rf_R'           ,...
                    'Ge'             ,...
					'Lz'             ,...
					'Zie'            ,...
					'Zfe'            ,...
                    'Tpr'            ,...
                    'se_flag'        ,...
                    'VOI_phi'        ,...
                    'VOI_theta'      ,...
                    'VOI_psi'        ,...
                    'VOI_X0'         ,...
                    'VOI_Y0'         ,...
                    'VOI_Z0'         ,...
                    'VOI_Lx'         ,...
                    'VOI_Ly'         ,...
                    'VOI_Lz'         ,...
                    'orient_flag'    ,...
                    'ovs_flag'       ,...
                    'OVS_FOV_x'      ,...
                    'OVS_FOV_y'      ,...
                    'OVS_FOV_z'      ,...
                    'OVS_delay'      ,...
                    'OVS_T'          ,...
                    'OVS_T_Kshift'   ,...
                    'OVStpwr'        ,...
                    'OVStpwrf'       ,...
                    'OVSsr'          ,...
                    'prfnum'         ,...
                    'tep'            ,...
                    'tDel'           ,...
                    'fix'};
                
str_params_list =  {'RFname'         ,...
                    'Gename'         ,...
                    'Ganame'         ,...
                    'OVSname'        ,...
                    'presig'};

fd = fopen(macro_fn,'w');
for idx = 1:length(real_params_list)
    fprintf(fd,'exists(''%s'',''parameter''):$dflag  if($dflag=1)  then  destroy(''%s'')  endif  create(''%s'',''real'')\n',...
               real_params_list{idx},real_params_list{idx},real_params_list{idx});
end;
for idx = 1:length(str_params_list)
    fprintf(fd,'exists(''%s'',''parameter''):$dflag  if($dflag=1)  then  destroy(''%s'')  endif  create(''%s'',''string'')\n',...
               str_params_list{idx},str_params_list{idx},str_params_list{idx});
end;

fprintf(fd,'\n'                                                       );
fprintf(fd,'Nsp=%f\n'            ,Nsp                                 );

% Chirp 90
fprintf(fd,'RFname=%s\n'         ,['''' chirp_name(1:end-3) '''']     );
fprintf(fd,'rfwdth=%f\n'         ,rfwdth                              );
fprintf(fd,'rfpwr=%1.0f\n'       ,rfpwr                               );
fprintf(fd,'rffpwr=%1.0f\n'      ,rffpwr                              );
fprintf(fd,'Gename=''%s''\n'     ,Gename                              );
fprintf(fd,'Ge=%f\n'             ,Ge                                  );
fprintf(fd,'Tp=%f\n'             ,Tp                                  );
fprintf(fd,'rf_R=%f\n'           ,rf_R                                );
fprintf(fd,'\n'                                                       );

% z-axis
fprintf(fd,'Zie=%f\n'            ,Zie                                 );
fprintf(fd,'Zfe=%f\n'            ,Zfe                                 );
fprintf(fd,'Lz=%f\n'             ,Lz                                  );

% Acquisition
fprintf(fd,'at=%f\n'             ,Ta                                  );
fprintf(fd,'Ga=%f\n'             ,Ga                                  );
fprintf(fd,'Ganame=''%s''\n'     ,Ganame                              );
fprintf(fd,'\n'                                                       );

% OVS
fprintf(fd,'VOI_phi=%3.3f\n'     ,0                                   );
fprintf(fd,'VOI_theta=%3.3f\n'   ,0                                   );
fprintf(fd,'VOI_psi=%3.3f\n'     ,0                                   );
fprintf(fd,'VOI_X0=%f\n'         ,0                                   );
fprintf(fd,'VOI_Y0=%f\n'         ,0                                   );
fprintf(fd,'VOI_Z0=%f\n'         ,0                                   );
fprintf(fd,'VOI_Lx=%f\n'         ,0                                   );
fprintf(fd,'VOI_Ly=%f\n'         ,0                                   );
fprintf(fd,'VOI_Lz=%f\n'         ,0                                   );
fprintf(fd,'ovs_flag=%f\n'       ,0                                   );
fprintf(fd,'OVS_FOV_x=%f\n'      ,6                                   );
fprintf(fd,'OVS_FOV_y=%f\n'      ,6                                   );
fprintf(fd,'OVS_FOV_z=%f\n'      ,8                                   );
fprintf(fd,'OVS_delay=%f\n'      ,10E-3                               );
fprintf(fd,'OVS_T=%f\n'          ,5E-3                                );
fprintf(fd,'OVS_T_Kshift=%f\n'   ,5E-3                                );
fprintf(fd,'OVStpwr=%f\n'        ,42                                  );
fprintf(fd,'OVStpwrf=%f\n'       ,4095                                );
fprintf(fd,'OVSsr=%f\n'          ,45E+3                               );
fprintf(fd,'OVSname=''%s''\n'    ,'ys_chirp90_T5_R10_SR45'            );
fprintf(fd,'\n'                                                       );

% Misc
fprintf(fd,'fix=%f\n'            ,1.0                                 );
fprintf(fd,'se_flag=%f\n'        ,se_flag                             );
fprintf(fd,'orient_flag=%f\n'    ,0                                   );
fprintf(fd,'ovs_flag=%f\n'       ,0                                   );
fprintf(fd,'prfnum=%f\n'         ,0                                   );
fprintf(fd,'phi=%3.1f\n'         ,0.0                                 );
fprintf(fd,'theta=%3.1f\n'       ,0.0                                 );
fprintf(fd,'psi=%3.1f\n'         ,0.0                                 );
fprintf(fd,'tep=%f\n'            ,tep                                 );
fprintf(fd,'tof=%f\n'            ,tof                                 );
fprintf(fd,'gain=%f\n'           ,gain                                );
fprintf(fd,'tDel=%f\n'           ,tDel                                );
fprintf(fd,'sw=%f\n'             ,sw                                  );
fprintf(fd,'presig=''l''\n'                                           );
fprintf(fd,'\n'                                                       );

fclose(fd);

% --------------------------------------------------------------------------
% FTP files
% --------------------------------------------------------------------------
declare_stage('FTPing files');
cur_pwd = pwd;

ftp_obj = ftp('frydmannmr300sun','imaging','varian1');
cd(root_dir);

if (exist(chirp_fn,'file') ~= 2)
    errordlg(sprintf('Cannot find 90-RF file (%s)',chirp_fn),'File Error');
else
    cd(ftp_obj,'/export/home/imaging/vnmrsys/shapelib');
    try
        disp(sprintf('FTPing %s\n',chirp_fn));
        mput(ftp_obj,chirp_fn);
    catch
        errordlg(sprintf('Failed to FTP 90-RF file (%s)',chirp_fn),'File Error');
        return;
    end;
end;

if (exist(macro_fn,'file') ~= 2)
    errordlg(sprintf('Cannot find macro file (%s)',macro_fn),'File Error');
else
    cd(ftp_obj,'/export/home/imaging/vnmrsys/maclib');
    try
        disp(sprintf('FTPing %s\n',macro_fn));
        mput(ftp_obj,macro_fn);
    catch
        errordlg(sprintf('Failed to FTP macro file (%s)',macro_fn),'File Error');
        return;
    end;
end;

close(ftp_obj);
cd(cur_pwd);

return;

