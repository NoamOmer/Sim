function [out_std_v,out_mean_std] = mean_std_around_extremum(mat,low_x,high_x,plot_f);
if (plot_f)
    figure;  subplot(1,3,1);  imagesc(mat);  title('Initial matrix');
    init_mat_sz = size(mat);
end;

mat = mat(low_x:high_x,:);
sz = size(mat);
N = sz(2);
x = (1:N);

mat = abs(mat);
if (plot_f)
    subplot(1,3,2);  imagesc(mat);
    axis([1,init_mat_sz(2),1-low_x,init_mat_sz(1)-low_x]);    title('Matrix abs value, after truncation');
end;

for idx = 1:sz(1)
    y = mat(idx,:);
    [yM,xM] = max(y);
    y = y.*(y>(0.05*yM));   % Thresh-hold

    tmp1 = (x.*y - xM*yM)./yM;
    tmp2 = (1/(N-1)) * sum(tmp1.^2);
    out_std_v(idx) = abs(xM - sqrt(tmp2));
end;
out_mean_std = mean(out_std_v);

if (plot_f)
    subplot(1,3,3); plot(out_std_v,(sz(1)-0):-1:(1-0)); title(sprintf('Matrix STD (mean=%4.4f)',out_mean_std));
    ylim([1-(init_mat_sz(1)-high_x),sz(1)+low_x]);
end;

