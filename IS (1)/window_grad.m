% Input
%   grad   : vector : [G/cm] : Time dependant gradient shape
%   t_grad : vector : [sec]  : Time axis
%   ESP    : scalar : [none] : Edge slope parameter
function [windowed_grad] = window_grad(grad,t_grad,ESP,DEBUG_FLAG);

tmp1 = length(grad);
tmp2 = [(-round(tmp1/2)+1):1:(+round(tmp1/2))];
GRDwin = exp(-(2.03*abs(tmp2./tmp1)).^(2*ESP));
GRDwin = GRDwin(1:tmp1);
windowed_grad = grad .* GRDwin;
if (DEBUG_FLAG >= 1)
	figure; hold;
	plot(t_grad, grad          ,'b.-');
	plot(t_grad, windowed_grad ,'g.-');
	title('Gradient(t)');
	xlabel('t_{grad} [sec]');
	ylabel('Gradient [G/cm]');
	legend({'pre-win','post-win'},'Location','Best');
	set_gca;
end;

return;

