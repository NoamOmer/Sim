
if (exist('context','var') == 1)  
	if (isempty(context))
		set_globals;
	else
		if (exist([context '.m'],'file') ~= 2)  
			error('Set_context: Error - missing context (%s.m)', context);
		else
			run(context);
		end;
	end;
else
	set_globals;
end;
