
% Finds the extremum point of a vector v
function [val,idx] = extrem(v)
[val,idx] = max(abs(v));
val = val * sign(v(idx));
return;

