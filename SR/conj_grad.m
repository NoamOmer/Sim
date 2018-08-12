% Conjugate Gradient Method Implementation with convergence check, finds
% the vector x for the solution of the linear equation Ax=b
%
% Implementation is based on: Algorithm 5.1, in page 73, in a PhD Thesis published by Stefan Kunis
% in 2006 and named "Nonequispaced FFT - Generalisation and Inversion"
%
% This version of the CG method assumes that the Ax=b problem is over-determined,
% hence, N, the size of the unknown vector x, is smaller than M, the size of the known b vector, i.e.,
%          f0  -->  x    size [N]       initial guess for the solution
%          S   -->  b    size [M]       known input
%          A   -->  A    size [MxN]     transformation matrix
% where,
%          N < M
%
% which implicitly allows the usage of non-square A matrices...

function [f,iteration_number,actual_accuracy, optimal_iterationNo] = conj_grad(f0,S,A,M,n_iter,epsilon)
set_globals;

actual_accuracy     = inf;
min_actual_accuracy = inf;
optimal_iterationNo = inf;

% Input data:
y  = S;                                                           % M

% Deviation of initial solution 'f0' from final solution 'f':
r0 = y - A*f0;                                                    % M = M - MxN * N

% Weight of the input data:
W  = ones(M,1);                                                   % M
% ...
p0 = ctranspose(A) * (W.*r0);                                     % N = NxM * (M.*M)

% ...
z0 = p0;                                                          % N

% Iteratively solve for f
for idx = 1:n_iter
	if (mod(idx,10) == 0)
	fprintf('idx = %1.0d  actual_accuracy = %+9.6g \t (epsilon=%-3.3g)\n',idx,actual_accuracy,epsilon);
	end;
	v0 = A*p0;                                                    % M = MxN * N
	a0 = (ctranspose(z0)*z0) / (ctranspose(v0)*(W.*v0));          % 1 = 1/1 = (N' * N) / (M' * (M.*M))
	f  = f0 + a0*p0;                                              % N
	r  = r0 - a0*v0;                                              % M

	new_actual_accuracy = (r'*r)/(S'*S);
	if (new_actual_accuracy < min_actual_accuracy)
		min_actual_accuracy = new_actual_accuracy;
		optimal_iterationNo = idx;
	end;
    actual_accuracy = new_actual_accuracy;
    if (actual_accuracy < epsilon)
        iteration_number = idx;
        return;            
    end
		
	z  = ctranspose(A) * (W.*r);                                  % N = NxM * (M.*M)
	b0 = (ctranspose(z)*z) / (ctranspose(z0)*z0);                 % 1 = 1/1 = (N' * N) / (N' * N)
	p  = z + b0*p0;                                               % N
	
	% Set the next iteration values
	f0 = f;
	r0 = r;
	p0 = p;
	z0 = z;
end;
iteration_number = idx;

return;

