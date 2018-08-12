function cc = set_cc(void)

cc = fix(clock);
dd = num2str(cc(3));
mm = num2str(cc(2));
yy = num2str(cc(1));  yy = yy(3:4);
cc = [dd,mm,yy,'_',num2str(cc(4)),num2str(cc(5)),num2str(cc(6))];

return;
