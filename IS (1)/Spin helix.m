clcl;
xlim = 0.8;
xres = 20e-3;
x = -xlim : xres : xlim ;

for loc = [0.01, 0.1, 0.2, 0.5, 0.99]
x0 = x(round(loc*length(x)));

a2 = 5;
c  = 0.5*exp(-5*abs(a2*(x-x0).^2));
c1 = 0.5 - c;
c2 = 0.5 - c;
c3 = 0.5 + c ;
% figure; plot(c1,'r-^'); hold on;
%         plot(c2,'b-+'); hold on;
%         plot(c3,'k.-'); legend({'c1','c2','c3'});
% return;
fh1 = figure; hold on;
for idx = 1:length(x)
	xi = x(idx);
	y  = a2*(xi-x0).^2;

	ye = exp(i*2*pi*y);
	a  = real(ye);
	b  = imag(ye);
	z  = 1e-5*xi;
	% figure; quiver(ones(1,length(a)),ones(1,length(a)),a,b);
	fh = quiver3(ones(1,length(a)),ones(1,length(a)),100*xi,a,b,z,25,'.');  hold on;
% 	figure(fh2);
	set(fh,'color',[c1(idx),c2(idx),c3(idx)],'LineWidth',2);
	% axis([0.9 1.1 0.9 1.1 -1.5 1.5]);
end;

x_HR = interp1(1:length(x),x,linspace(1,length(x),20*length(x)));
y  = a2*(x_HR -x0).^2;
ye = exp(i*2*pi*y);
a  = real(ye);
b  = imag(ye);
z  = 100*x_HR;
ax = plot3(26*a,26*b,z,'-','LineWidth',1); hold on;
set(ax,'color',[0.5 0.5 0.5],'MarkerSize',4);

grid off;
axis image;
set(gca,'xtick',[],'ytick',[],'ztick',[]);                    % remove the axes ticks
set(gca,'view',[-17,12]);                                     % set view angle
set(gca,'XColor',[1,1,1],'YColor',[1,1,1],'ZColor',[1,1,1]);  % blank out the axes frame

end; % loc