% Rotation around y axis in angle phi

function Ry   = yrot(phi)

Ry = [cos(phi) 0 sin(phi); 0 1 0; -sin(phi) 0 cos(phi)];

return;

