% Elliptical approximation using bisection for c search (broken)

function [Xf_H, Xf_h, c_opt] = calculate_Xf(A, B, K, P, F, E, e, xr, ur)

    % Initial upper/lower bounds
    c_low = 0;
    c_high = 1e5;
    tol = 1e-3;
    max_iter = 50;

    F_tilde = F + E*K;
    e_bar = e - F*xr - E*ur;
    
    n = size(A, 1);
    c_opt = 0;

    for iter = 1:max_iter
        c_mid = (c_low + c_high) / 2;
        disp(['Trying c = ', num2str(c_mid)])

        polyhedron = compute_scaled_ellipsoid_polyhedron(P, c_mid);
        vertices = polyhedron.V;

        feasible = true;
        for i = 1:size(vertices, 1)
            xbar = vertices(i, :)';
            xbar_next = (A + B*K) * xbar;
            ubar = K * xbar;

            % Check if successor is still in the set
            if ~all(polyhedron.A * (xbar_next) <= polyhedron.b + ones(size(polyhedron.b))*1e-8)
                disp("Failed invariance")
                % % disp(x)
                % disp(max(polyhedron.A * (x_plus) - polyhedron.b))
                % % disp(polyhedron.A * (x_plus) - polyhedron.b)
                feasible = false;
                break
            end
            % Check if constraints are satisfied
            if ~all(F*(xbar+xr) + E * (ubar+ur) <= e + ones(size(e))*1e-8)
                disp("Failed constraints")
                % disp(F_tilde*(xbar+xr) <= e)
                % disp(F*(x) + E*(u))
                feasible = false;
                break
            end
        end

        if feasible
            % increase c
            c_low = c_mid;
            c_opt = c_mid;
            best_poly = polyhedron;
        else
            % decrease c
            c_high = c_mid;
        end

        % Converged
        if abs(c_high - c_low) < tol
            break
        end
    end

    if c_opt > 0
        best_poly = best_poly + xr;
        Xf_H = best_poly.A;
        Xf_h = best_poly.b;
    else
        warning('No feasible terminal set found.');
        Xf_H = [];
        Xf_h = [];
    end
end
% 
% 
function polyhedron = compute_scaled_ellipsoid_polyhedron(P, c)

    % Step 1: Eigen decomposition of P
    [V, ~] = eig(P);  % V contains eigenvectors (orthonormal)

    % Step 2: Use directions along each eigenvector (and their negatives)
    directions = [V, -V];  % Each column is a direction vector

    random_directions = randn(size(P,1), 50);
    random_directions = random_directions ./ vecnorm(random_directions);
    directions = [directions random_directions];

    N = size(directions, 2);

    % Step 3: Compute support function for each direction
    A_poly = zeros(N, size(P,1));
    b_poly = zeros(N, 1);

    for i = 1:N
        d = directions(:, i);

        % Support function: h(d) = sqrt(c * d' * P^{-1} * d)
        h_i = sqrt(c * (d' * (P \ d)));%sqrt(c) * sqrt(d' * inv(P) * d) + d'*xr;  % Equivalent to: d' * inv(P) * d

        A_poly(i, :) = d';
        b_poly(i) = h_i;
    end
    polyhedron = Polyhedron('A', A_poly, 'b', b_poly);
end

% % Maximal control invariant set
% function [Xf_H, Xf_h] = calculate_Xfb(A, B, K, F, E, e, xr, ur)
% system = LTISystem('A', A, 'B', B);
% 
% % Translate constraints to xr and ur
% F_tilde = (F + E*K);  % Translate constraints to only x by using u = K*x
% e_bar = e + F*xr - E*ur;
% 
% % Make constraints
% 
% constraints_orig = Polyhedron('A', F_tilde, 'b', e);
% system.setDomain('x', constraints_orig);
% 
% % Control invariant set
% global MPTOPTIONS
% MPTOPTIONS = mptopt('verbose', 2);
% 
% tic;
% Xf = system.invariantSet();
% toc;
% Xf_H = Xf.A;
% Xf_h = Xf.b;
% 
% end
% 
% 
% 
% % Comparing my elliptical approximation with maximal invariant set example from slides (Lec. 4 10-11)
% 
% A = [1 1;
%     0 1];
% B = [0;
%     1];
% Q = eye(2);
% R = 3;
% 
% xr = [0.0;
%      0.3];
% ur = 0.0;
% [P,K,L] = idare(A,B,Q,R); K=-K;
% 
% F = [ 1  0;   %  x1 <= 3/4
%      -1  0;   % -x1 <= 3/4
%       0  1;   %  x2 <= 3/4
%       0 -1;   % -x2 <= 3/4
%       0  0;   %  no x term in u constraints
%       0  0];  %  no x term in u constraints
% E = [ 0;      %  no u term in x constraints
%       0;
%       0;
%       0;
%       1;      %  u <= 1/2
%      -1];     % -u <= 1/2
% e = [3/4;  % x1 upper bound
%      3/4;  % x1 lower bound 
%      3/4;  % x2 upper bound
%      0;  % x2 lower bound
%      1/2;  % u upper bound
%      1/2]; % u lower bound
% 
% [Xf_H, Xf_h] = calculate_Xfa(A, B, K, P, F, E, e, xr, ur);
% my_inv_set = Polyhedron('A', Xf_H, 'b', Xf_h);
% [Xf_H, Xf_h] = calculate_Xfb(A, B, K, F, E, e, xr, ur);
% my_inv_set2 = Polyhedron('A', Xf_H, 'b', Xf_h);
% 
% F_tilde = F + E*K;
% constraints = Polyhedron('A', F_tilde, 'b', e);
% figure(1);
% plot(constraints, 'alpha', 0.1, 'color', 'red')
% hold on
% plot(my_inv_set2, 'color', 'green', 'alpha', 0.2)
% plot(my_inv_set, 'color', 'blue', 'alpha', 0.3);
% hold off
% legend('Constraints', 'Terminal set (maximal)', 'Terminal set (ellipse)')
% % hold off
% % axis equal
% 
% %% My own example from our system
% 
% [A,B,~] = ltvStockModel(50);
% Q = 1*eye(5);
% R = 2*eye(8);
% [P,K,L] = idare(A,B,Q,R);K=-K;
% xr = [1, 1, 1, 1, 200]';
% ur = ones(8,1);
% 
% [xfH, xfh] = calculate_Xfa(A, B, K, P, F, E, e, xr, ur);
