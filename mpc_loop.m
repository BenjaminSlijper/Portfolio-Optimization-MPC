function [u] = mpc_loop(x0,x_ref, Q, R, F, E, e, k, N)

alpha = 0.0001; 
beta = 0.0001;
n = 4;

u_val = sdpvar(2*n, N);  % nu is the number of control inputs, N is the prediction horizon

cost = 0;
constraints = [];

for i = 1:N
    [A,B,C,D,E_,dx0,x0_,u0,y0,Delays] = ltvStockModel(i+k);
    cost = cost + (x0-x_ref)'*Q*(x0-x_ref) + u_val(:,i)'*R*u_val(:,i);
    constraints = [constraints, F*x0 + E*u_val(:,i) <= e];
    x0 = A*x0 + B*u_val(:,i);
end
P = dare(A,B,Q,R);

cost = cost + (x0-x_ref)'*P*(x0-x_ref);

options = sdpsettings('verbose', 0);
optimize(constraints,cost);

u = value(u_val);