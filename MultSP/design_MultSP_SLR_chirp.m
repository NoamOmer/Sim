
function [chirp_rot,rfpwr,rffpwr] = design_MultSP_SLR_chirp(GePE,y_axis,dte,Nsp,alpha2,Tp,rfwdth,...
                                                            tpwr90,pw90,RF_ESP_SLR)
set_globals;

OmegaE      = gammaHz*GePE*y_axis;
OmegaE_axis = OmegaE;
y_of_OmegaE = interp1(OmegaE,y_axis,OmegaE);
swHz        = 1/dte;
isPI        = 0;
if (Nsp > 1)
	n = round(linspace(0,length(y_axis),Nsp+1));
	Phi_e_sq_rad = [];
	for sp_idx = 1:Nsp
		ax = y_axis(n(sp_idx)+1 : n(sp_idx+1));                               % Spatial range of current SP
		mid_val = mean([y_axis(n(sp_idx)+1), y_axis(n(sp_idx+1))]);           % Middle value of the range
		Phi_e_sq_rad(end+1:n(sp_idx+1)) = 2*pi*alpha2*((ax - mid_val) .^ 2);  % Frequency response
	end;
else
	Phi_e_sq_rad = 2*pi*alpha2*(y_axis.^2);
end;

if (DEBUG_FLAG >= 5)
figure;
subplot(3,1,1); plot(y_axis     ,OmegaE*1E-3 ,'.-'); title('OmegaE(y)  '); xlabel('y-axis [cm]'   ); ylabel('Omega(y) [kHz]'       ); grid;
subplot(3,1,2); plot(OmegaE*1E-3,y_of_OmegaE ,'.-'); title('y(OmegaE)  '); xlabel('Omega(y) [kHz]'); ylabel('y-axis [cm]'          ); grid;
subplot(3,1,3); plot(y_axis     ,Phi_e_sq_rad,'.-'); title('Phi e sq(y)'); xlabel('y-axis [cm]'   ); ylabel('Phi_e_sq_rad(y) [rad]'); grid;
end;
[chirp_rot,rfpwr,rffpwr] = convert_pulse_to_SLR(y_axis,Tp,rfwdth,Phi_e_sq_rad,OmegaE,OmegaE_axis,...
	                                            y_of_OmegaE,swHz,isPI,tpwr90,pw90,RF_ESP_SLR,...
											    DEBUG_FLAG,SIMULATION);

return;
