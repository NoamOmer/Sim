function OS_compatible_path = nbePath(p,slash_input)

if (isempty(p) || (~isstr(p)))
	return;
end;

loc1 = findstr(p,'/');
loc2 = findstr(p,'\');

if (exist('slash_input','var'))
	p(loc1) = slash_input;
	p(loc2) = slash_input;
else
	p(loc1) = filesep;
	p(loc2) = filesep;
end;

OS_compatible_path = p;

return
