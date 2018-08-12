% Rotation around z axis in angle phi - counter-clockwise rotation

function Rz = zrot(phi)

Rz = [cos(phi) -sin(phi) 0;sin(phi) cos(phi) 0; 0 0 1];

return;
