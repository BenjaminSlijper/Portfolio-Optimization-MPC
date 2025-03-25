%Define Model Parameters

alpha = 0.0001; 
beta = 0.0001;
n = 4;
N = 10;
Nsim = 10;

[F,E,e]=constraintgen(n, alpha, beta); %Generate Model Constraints

Q = ones(5, 5);  % Example state cost weighting matrix
R = eye(2*n);  % Example input cost weighting matrix

x = zeros(5,Nsim);
x(:,1) = [0;0;0;0;100];
x_ref = [1000;1000;1000;1000;1000];
y = zeros(Nsim,1);

for k = 1:Nsim
    u = mpc_loop(x(:,k), x_ref, Q, R, F, E, e, k, N);
    [A,B,C,D,E_,dx0,x0_,u0,y0,Delays] = ltvStockModel(k);
    x(:,k+1) = A*x(:,k)+B*u(:,1);
    y(k,:) = C*x(:,k);
end

%Plotting

figure;
plot(x(1:4,:)');
title('Value of Each Stock Owned');
xlabel('Time (days)');
ylabel('Stock Value (USD)');
legend('Stock 1', 'Stock 2', 'Stock 3', 'Stock 4'); % Adjust if necessary
grid on;


figure;
plot(y);
title('Total Portfolio Value');
xlabel('Time (days)');
ylabel('Value (USD)');
grid on;

