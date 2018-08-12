function [m_rot] = nbe_rot_mat(m,deg)

[mx,my]   = size(m);
 m1       = imrotate(m,deg,'nearest','loose');
[m1x,m1y] = size(m1);
start_x   = round(m1x/2 - mx/2);
start_y   = round(m1y/2 - my/2);
m_rot     = m1( (start_x : start_x + mx - 1) , (start_y : start_y + my - 1) );

return;

