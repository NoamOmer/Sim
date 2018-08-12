% This function produces a vector rotation matrix. Given a vector v
% the matrix R will rotate the vector around the arbitry unit vector
% in angle phi using the LEFT HAND RULE
% R*(v') = rotated vector
% This function does not rotate the coordinate system but the vector.

% Positive angle + positive unit vector -->         clockwise rotation 
% Positive angle + negative unit vector --> counter-clockwise rotation
% Negative angle + positive unit vector --> counter-clockwise rotation        
% Negative angle + negative unit vector -->         clockwise rotation 

function R = rot(unit_vec,phi)

phi = -phi;
x = unit_vec(1);
y = unit_vec(2);
z = unit_vec(3);
if ( (abs(sqrt(sum(unit_vec.^2)) - 1) > 1E-13) || isnan(x) || isnan(y) || isnan(z) )
	error('rot: x,y,z do not represent a unit vector (x=%d, y=%d, z=%d). Exiting',x,y,z);
end;

t = 1 - cos(phi);
s = sin(phi);
c = cos(phi);

R = [t*(x^2) + c,  t*x*y - s*z,  t*x*z + s*y; ...
     t*x*y + s*z,  t*(y^2) + c,  t*y*z - s*x; ...
     t*x*z - s*y,  t*y*z + s*x,  t*(z^2) + c];

return;

