%Definition of system dimension
dim.nx=5;     %state dimension
dim.ny=1;     %output dimension
dim.nu=8;     %input dimension
dim.N=20;      %horizon

Q = eye(dim.nx, dim.nx);
R = eye(dim.nu, dim.nu);
P = 10 * Q;  % make better P


Nsim = 1200;
x0 = [0, 0, 0, 0, 100]';
x = zeros(dim.nx,Nsim);
x(:,1)=x0;
y = zeros(1, Nsim);
y(1) = sum(x0);

u_rec=zeros(dim.nu,Nsim+1);

for k=1:Nsim
    % replace with optimal state space solution
    xr = [100, 100, 100, 100, 100]';
    ur = [0, 0, 0, 0, 0, 0, 0, 0]';

    [T,S] = predmodgen_ltv(@ltvStockModel, k, dim);
    [H,h,const] = costgen(T, S, Q, R, P, dim, x(:, k), xr, ur);
    
    %Solve optimization problem    
                                  %define optimization variable
    % Constraint=[];                                                 %define constraints
    % uostar = sdpvar(dim.nu*dim.N,1);
    % Objective = 0.5*uostar'*H*uostar+h'*uostar;    %define cost function
    % optimize(Constraint,Objective);                                %solve the problem
    % uostar=value(uostar);

    uostar = quadprog(H, h);

    % Select the first input only
    u_rec(:,k)=uostar(1:dim.nu);
    
    [A,B,C,~] = ltvStockModel(k);
    % Compute the state/output evolution
    x(:,k+1)=A*x(:,k) + B*u_rec(:,k);
    y(k+1)=C*x(:,k+1);
end

%%

for i=1:5
    figure(i);
    plot(0:Nsim, x(i, :))
    title((sprintf("Stock %i", i)))
    hold on
    yline(100, "Color", 'black')
    hold off
end
figure(6);
plot(0:Nsim, y)
title("Total portfolio")
