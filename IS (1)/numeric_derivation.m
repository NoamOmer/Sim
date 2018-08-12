% Numeric derivation of the form dy/dx

function [dy_dx] = numeric_derivation(y,x);

sz = size(y);
if ((sz(1) == 1) || (sz(2) == 1))
	% Calculate dO/dte
	dy_dx = diff(y)./(x(2)-x(1));     % decreases the vector length by 1
	dy_dx(end+1) = dy_dx(end);        % resize to length of 'te'
else
	for idx = 1:sz(1)
		% Hey! this is not efficient, it can be done using diff(MAT,1,1) or diff(MAT,1,2)
		dy_dx(idx,1:sz(2)) = numeric_derivation(y(idx,:),x);
	end;
end;

return;

