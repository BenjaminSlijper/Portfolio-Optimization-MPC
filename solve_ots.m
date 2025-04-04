function [xr,ur] = solve_ots(A, B, C, R, sigma, y_ref, F, E, e, dim)
    lambda = 100;
    Q = eye(dim.nu);

    x_ref = sdpvar(dim.nx, 1);
    u_ref = sdpvar(dim.nu, 1);

    Objective = -R'*x_ref + lambda * x_ref' * sigma * x_ref + u_ref' * Q * u_ref;
    
    steady_state_constraint = [x_ref == A*x_ref + B*u_ref];
    ref_constraint = [y_ref == C*x_ref];
    state_constraints = [F * x_ref + E * u_ref <= e];
    constraints = [steady_state_constraint, ref_constraint, state_constraints];
    % constraints = [A_eq * [x_ref; u_ref] == b_eq];
    % constraints = [];%[sum(x_ref) <= y_ref];
    options = sdpsettings('solver', 'quadprog', 'verbose', 0);
    optimize(constraints, Objective, options);

    xr = value(x_ref);
    ur = value(u_ref);
end

% function [xr, ur] = solve_ots(A, B, C, R, sigma, y_ref, dim)
%     % Parameters
%     lambda = 1;  % Risk aversion parameter
%     gamma = 0.1; % Trading cost parameter
%     Q = eye(dim.nu);
% 
%     % Decision variables
%     x_ref = sdpvar(dim.nx, 1);
%     u_ref = sdpvar(dim.nu, 1);
% 
%     % Objective: maximize return, minimize risk and control effort
%     Objective = -R'*x_ref + lambda*x_ref'*sigma*x_ref + gamma*u_ref'*Q*u_ref;
% 
%     % Constraints
%     constraints = [];
%     % Steady-state constraint
%     constraints = [constraints, x_ref == A*x_ref + B*u_ref];
%     % Output tracking constraint
%     % constraints = [constraints, C*x_ref == y_ref];
%     % Optional: Add position constraints if needed
%     constraints = [constraints, x_ref >= 0];
% 
%     % Solve
%     options = sdpsettings('verbose', 1);
%     solution = optimize(constraints, Objective, options);
% 
%     % Return optimal values
%     xr = value(x_ref);
%     ur = value(u_ref);
% 
%     % Debug
%     if solution.problem ~= 0
%         disp(['Optimization failed: ', yalmiperror(solution.problem)]);
%     else
%         disp(['Target return: ', num2str(y_ref)]);
%         disp(['Achieved return: ', num2str(C*xr)]);
%         disp(['Portfolio risk: ', num2str(xr'*sigma*xr)]);
%     end
% end