yalmip('clear')

%% Variables

%Definition of system dimension
dim.nx=5;     %state dimension
dim.ny=1;     %output dimension
dim.nu=8;     %input dimension
dim.N=20;      %horizon

Nsim = 1000;  %simulation time

x0 = [0, 0, 0, 0, 100]';  %initial condition
y_ref = 200;  %TODO find optimal y_ref in the MPC loop

Q = eye(dim.nx, dim.nx);
R = eye(dim.nu, dim.nu);

alpha = 0.0001;
beta = 0.0001;

% --- Compute covariance matrix and return matrix
total_return_matrix = zeros(1200, dim.nx);
for k=1:1200
    [~, return_matrix] = ltvStockModel(k);
    total_return_matrix(k, :) = diag(return_matrix);
end
sigma = cov(total_return_matrix);
avg_returns = mean(total_return_matrix, 1)';

%% MPC

% Setting up variables
x = zeros(dim.nx,Nsim);  %state trajectory
x(:,1)=x0;
y = zeros(1, Nsim);  %output trajectory
y(1) = sum(x0);

u_rec=zeros(dim.nu,Nsim+1);  %optimal inputs over time

for k=1:Nsim
    disp(k)
    
    % Setting up matrices
    [T,S] = predmodgen_ltv(@ltvStockModel, k, dim);
    [A_ineq,b_ineq, F,E,e] = constraintgen(4, alpha, beta, T, S, dim, x(:, k));

    [A_N,B_N,C_N] = ltvStockModel(k+dim.N);  %state matrices at end of prediction horizon
    [xr,ur] = solve_ots(A_N, B_N, C_N, avg_returns, sigma, y_ref, F, E, e, dim);

    [P,~] = idare(A_N, B_N, Q, R);
    [H,h,const] = costgen(T, S, Q, R, P, dim, x(:, k), xr, ur);

    %Solve optimization problem    
    warning off all
    options = optimoptions('quadprog', 'Display', 'off');
    uostar = quadprog(H, h, A_ineq, b_ineq, [], [], [], [], [], options);
    warning on all

    % Select the first input only
    u_rec(:,k) = uostar(1:dim.nu);

    % Compute the state/output evolution
    [A,B,C,~] = ltvStockModel(k);
    x(:,k+1) = A*x(:,k) + B*u_rec(:,k);
    y(k+1) = C*x(:,k+1);
end

%% Plotting

for i=1:5
    figure(i);
    plot(0:Nsim, x(i, :))
    title((sprintf("Stock %i", i)))
end
figure(6);
plot(0:Nsim, y)
hold on
yline(y_ref, "Color", 'black')
hold off
title("Total portfolio")
