
function [f] = increase_1Dvector_res(f0,S,A,M,x,n_iter)
% set_globals; - not needed and SLOWS things down!

% Input data:
% % % % % y  = S;        commented for speed                                                   % M

% Deviation of initial solution 'f0' from final solution 'f':
% % % % % r0 = y - A*f0; commented for speed                                                   % M = M - MxN * N
r0 = S - A*f0;                                                    % M = M - MxN * N

% Weight of the input data:
W  = ones(M,1);                                                   % M
% ...
p0 = ctranspose(A) * (W.*r0);                                     % N = NxM * (M.*M)

% ...
z0 = p0;                                                          % N

% Iteratively solve for f
% % % % % for idx = 1:n_iter
	v0 = A*p0;                                                    % M = MxN * N
	a0 = (ctranspose(z0)*z0) / (ctranspose(v0)*(W.*v0));          % 1 = 1/1 = (N' * N) / (M' * (M.*M))
	f  = f0 + a0*p0;                                              % N

% 	-------------------------------------
% 	comment the following block for speed
%                  (and the forloop)
% 	-------------------------------------
% % % % % 	r  = r0 - a0*v0;                                              % M
% % % % % 	z  = ctranspose(A) * (W.*r);                                  % N = NxM * (M.*M)
% % % % % 	b0 = (ctranspose(z)*z) / (ctranspose(z0)*z0);                 % 1 = 1/1 = (N' * N) / (N' * N)
% % % % % 	p  = z + b0*p0;                                               % N
% % % % % 	
% % % % % 	% Set the next iteration values
% % % % % 	f0 = f;
% % % % % 	r0 = r;
% % % % % 	p0 = p;
% % % % % 	z0 = z;
% % % % % 
% % % % % % 	if (DEBUG_FLAG >= 5)
% % % % % % 	figure; plot(linspace(x(1),x(end),M), abs(S) / max(abs(S)),'k.-'); hold on;
% % % % % % 	        plot(x                      , abs(f) / max(abs(f)),'r.-'); legend({'Original S','SR S'});
% % % % % % 	        title(sprintf('SR profile (iter = %1.0f)',idx));
% % % % % % 	end;
% % % % % end;

return;

