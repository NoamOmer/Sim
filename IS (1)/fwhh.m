% Finds the full width at half height of an input vector

function [fwhh_val] = fwhh(x,y,interp_factor,DEBUG_FLAG)
xi = x(1) : ((x(2)-x(1))/interp_factor) : x(end);
yi = interp1(x,y,xi);

max_yi = max(yi);
[dummy fwhh_idx] = sort(abs(yi - max_yi/2));
ctr  = 1;
idx1 = fwhh_idx(ctr);
idx2 = idx1;
while (abs(idx1-idx2) < 5)
    ctr = ctr + 1;
    idx2 = fwhh_idx(ctr);
end;
fwhh_val = abs(xi(idx1) - xi(idx2));

% Plot
if (DEBUG_FLAG >= 3)
    fh = figure; hold on;
    plot(x,y,'ko',xi,yi,'m.-');

    left  = min(idx1,idx2);
    right = max(idx1,idx2);
    line([xi(left) xi(right)],[yi(left) yi(right)]);
    legend({'Pre-interp','Post-interp'});
end;
