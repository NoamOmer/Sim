function changeDisplayInterval
u=gcbo;
Sli=get(u,'userdata');
a=inputdlg('');
if isempty(a),return;end
a=str2num(a{1});
U=get(Sli,'userdata');
U.displayInterval=a;
U.minStep=a/10;
U.maxStep=a-1;
set(Sli,'userdata',U);
u=get(gca,'userdata');
len=u{6}-U.SF*a;
set(Sli,'max',len);
set(Sli,'sliderstep',[a/10 a*.9]*U.SF/len);