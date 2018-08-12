function plot_square(x,y,fmt,lnwidth)
plot([x(1) x(1)],[y(1) y(2)],fmt); a = get(gca,'children'); set(a(1),'LineWidth',lnwidth);
plot([x(2) x(2)],[y(1) y(2)],fmt); a = get(gca,'children'); set(a(1),'LineWidth',lnwidth);
plot([x(1) x(2)],[y(1) y(1)],fmt); a = get(gca,'children'); set(a(1),'LineWidth',lnwidth);
plot([x(1) x(2)],[y(2) y(2)],fmt); a = get(gca,'children'); set(a(1),'LineWidth',lnwidth);
return;
