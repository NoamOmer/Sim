% STD according to FWHH (Full Width @ Half Height)
% Function calculates the value around an input matrix extremum
% meaning that it will use the absolute value of the matrix.
% (not recommended for complex matrices).
function [out_std_v,out_mean_std] = mean_fwhh_std(mat,low_x,high_x,plot_f);
if (plot_f)
    figure;  subplot(1,3,1);  imagesc(mat);  title('Initial matrix');
    init_mat_sz = size(mat);
end;

% Resize matrix
mat = mat(low_x:high_x,:);
sz  = size(mat);
N   = sz(2);
x   = (1:N);

mat = abs(mat);
if (plot_f)
    subplot(1,3,2);  imagesc(mat);
    axis([1,init_mat_sz(2),1-low_x,init_mat_sz(1)-low_x]);    title('Matrix abs value, after truncation');
end;

for idx = 1:sz(1)
    INTERP_FACTOR = 10;
	y = mat(idx,:);
    x = 1:length(y);                                   % Create x-axis (just indices)
    xi = x(1) : ((x(2)-x(1))/INTERP_FACTOR) : x(end);  % Interpulate x-axis
    yi = interp1(x,y,xi);                              % Interpulate y-axis
    yi = smooth(yi);                                   % Smooth y-axis
    
    [yMax,xMax] = max(yi);                             % Find the maximal y-value
    xright = find((yi(xMax:end) - yMax/2) < 0) - 1;    % Move to the right and stop at half height
    xright = xMax + xright(1);                         % Take the first index that passes the half height
    xleft  = find((yi(xMax:-1:1) - yMax/2) < 0);
    xleft  = xMax - xleft(1);
    out_std_v(idx) = round((xright - xleft)/INTERP_FACTOR);
%    [dummy,xleft]  = min(abs(y(1:xMax)   - yMax/2));
%    [dummy,xright] = min(abs(y(xMax:end) - yMax/2));
%    out_std_v(idx) = xMax + xright - xleft;
end;
out_mean_std = mean(out_std_v);

if (plot_f)
    subplot(1,3,3); plot(out_std_v,sz(1):-1:1); title(sprintf('Matrix STD3 (mean=%4.4f)',out_mean_std));
    xlabel('FWHH [pixels]');
    ylim([1-(init_mat_sz(1)-high_x),sz(1)+low_x]);
end;

