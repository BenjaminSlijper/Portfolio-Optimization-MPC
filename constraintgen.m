function [F,E,e]=constraintgen(n, alpha, beta)

F = [zeros(2*n, n+1) ; zeros(1,n) -1 ; -eye(n) zeros(n,1); zeros(1,n) -1];
E = [-eye(2*n) ; (1+alpha)*ones(1,n) -(1+beta)*ones(1,n) ; -eye(n) eye(n); zeros(1,2*n)];
e = zeros(3*n+2,1);
