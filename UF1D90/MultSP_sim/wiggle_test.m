clcl;
load wiggle_test;
p = Phi_M;
p1 = p(250:2150);
f1 = exp(i*p1);
figure;
subplot(2,1,1); plot(abs(cumsum(f1)),'.-'); title('f1');
subplot(2,1,2); plot(           p1  ,'.-'); title('p1');

inv_p1 = p1(end:-1:1);
inv_f1 = exp(i*inv_p1);
figure;
subplot(2,1,1); plot(abs(cumsum(inv_f1)),'.-'); title('inv f1');
subplot(2,1,2); plot(           inv_p1  ,'.-'); title('inv p1');

p1_left  = p1(1:983);
p1_right = p1(984:end);
f1_left  = exp(i*p1_left);
f1_right = exp(i*p1_right);
figure;
subplot(2,1,1); plot(abs(cumsum(f1_left )),'.-'); title('f1 left ');
subplot(2,1,2); plot(           p1_left   ,'.-'); title('p1 left ');
figure;
subplot(2,1,1); plot(abs(cumsum(f1_right)),'.-'); title('f1 right');
subplot(2,1,2); plot(           p1_right  ,'.-'); title('p1 right');

inv_p_left  = p1(983:-1:1);
inv_p_right = p1(end:-1:984);
inv_f_left  = exp(i*inv_p_left);
inv_f_right = exp(i*inv_p_right);
figure;
subplot(2,1,1); plot(abs(cumsum(inv_f_left )),'.-'); title('inv f left ');
subplot(2,1,2); plot(           inv_p_left   ,'.-'); title('inv p left ');
figure;
subplot(2,1,1); plot(abs(cumsum(inv_f_right)),'.-'); title('inv f right');
subplot(2,1,2); plot(           inv_p_right  ,'.-'); title('inv p right');
