
% Anti clockwise rotation of a point (x,y) around the origin (0,0)
% at an angle phi [rad]

function [rotated_x, rotated_y] = rot_xy_Point(x,y,phi)

rotated_x = x * cos(phi) - y * sin(phi);
rotated_y = x * sin(phi) + y * cos(phi);

return;

