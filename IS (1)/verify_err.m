function verify_err(s,err_lvl)
% uiwait(msgbox(s));
rc = questdlg([s '...  continue?'],'Verification error','Yes','No','Yes');
if strcmp(rc,'No')
	error(s);
else
	disp('***');
	disp(s);
	disp('***');
end;
return;

