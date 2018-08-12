function [cdata,x,y] = get_fig_data(is_subplot_f, plot_f)

ch=get(gca,'children');
if (exist('is_subplot_f','var') && ~isempty(is_subplot_f) && (is_subplot_f~=0))
	ch = get(gcf,'children');
	ch = ch(1);
	ch = get(ch,'children');
	ch = ch(is_subplot_f);
	x  = get(ch,'Xdata');
	y  = get(ch,'Ydata');
	cdata=y;
	figure;plot(x,y,'.-');
else
	cdata=get(ch,'CData');
end;

if exist('plot_f','var') && plot_f
	figure; imagesc(cdata);
end;
return;