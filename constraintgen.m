function [A,b]=constraintgen(n, alpha, beta, T, S, dim, x0)
% what is this n? is it n_inputs/2?
F = [zeros(2*n, n+1) ; zeros(1,n) -1 ; -eye(n) zeros(n,1); zeros(1,n) -1];
E = [-eye(2*n) ; (1+alpha)*ones(1,n) -(1+beta)*ones(1,n) ; -eye(n) eye(n); zeros(1,2*n)];
e = zeros(3*n+2,1);

F_list = repmat({F}, dim.N, 1);
E_list = repmat({E}, dim.N, 1);

F_big = blkdiag(F_list{:});
E_big = blkdiag(E_list{:});
e_big = repmat(e, dim.N, 1);

A = F_big * S + E_big;
b = e_big - F_big * T * x0;
