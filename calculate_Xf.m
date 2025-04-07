% function [Xf_H, Xf_h] = calculate_Xfa(A, B, K, P, F, E, e)
% c = 10; 
% feasible = false;
% 
% while ~feasible
%     c = c / 1.01;
%     disp(c)
%     polyhedron = compute_scaled_ellipsoid_polyhedron(P, c);
% 
%     % mapped_poly = polyhedron.affineMap(A + B*K);
%     % control_invariant = polyhedron.contains(mapped_poly);
% 
% 
% 
%     vertices = polyhedron.V;
%     A_poly = polyhedron.A;
%     b_poly = polyhedron.b;
% 
%     control_invariant = true;
%     for i=1:length(polyhedron.V)
%         vertex = vertices(i, :)';
%         u = K*vertex;
%         x_plus = (A + B*K)*vertex;
% 
%         if ~all(A_poly * x_plus <= b_poly)
%             control_invariant = false;
%             break
%         end
%         if ~all(F*x + E*u <= e)
%             control_invariant = false;
%             break
%         end
%     end
% 
%     if control_invariant
%         feasible = true;
%         Xf_H = polyhedron.A;
%         Xf_h = polyhedron.b;
%     end
% end
% end

function [Xf_H, Xf_h, c_opt] = calculate_Xfa(A, B, K, P, F, E, e, x_ref, u_ref)

    % Initial upper/lower bounds
    c_low = 0;
    c_high = 1e4;
    tol = 1e-3;
    max_iter = 50;
    
    n = size(A, 1);
    c_opt = 0;

    for iter = 1:max_iter
        c_mid = (c_low + c_high) / 2;
        disp(['Trying c = ', num2str(c_mid)])

        polyhedron = compute_scaled_ellipsoid_polyhedron(P, c_mid);
        vertices = polyhedron.V;

        feasible = true;
        for i = 1:size(vertices, 1)
            x = vertices(i, :)';
            u = K * x;
            x_plus = A*x + B*u;

            % Check if successor is still in the set
            if ~all(polyhedron.A * (x_plus) <= polyhedron.b)
                % disp("Failed invariance")
                % disp(x)
                % disp(polyhedron.A * (x_plus) <= polyhedron.b)
                feasible = false;
                break
            end
            % Check if constraints are satisfied
            if ~all(F*(x+x_ref) + E*(u+u_ref) <= e)
                % disp("Failed constraints")
                % disp(x)
                % disp(F*(x+x_ref) + E*(u+u_ref))
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
        Xf_H = best_poly.A;
        Xf_h = best_poly.b;
    else
        warning('No feasible terminal set found.');
        Xf_H = [];
        Xf_h = [];
    end
end


function polyhedron = compute_scaled_ellipsoid_polyhedron(P, c)

    % Step 1: Eigen decomposition of P
    [V, ~] = eig(P);  % V contains eigenvectors (orthonormal)

    % Step 2: Use directions along each eigenvector (and their negatives)
    directions = [V, -V];  % Each column is a direction vector
    
    % random_directions = randn(2, 100);
    % random_directions = random_directions ./ vecnorm(random_directions);
    % directions = [directions random_directions];

    N = size(directions, 2);

    % Step 3: Compute support function for each direction
    A_poly = zeros(N, size(P,1));
    b_poly = zeros(N, 1);

    for i = 1:N
        d = directions(:, i);

        % Support function: h(d) = sqrt(c * d' * P^{-1} * d)
        h_i = sqrt(c * (d' * (P \ d)));  % Equivalent to: d' * inv(P) * d

        A_poly(i, :) = d';
        b_poly(i) = h_i;
    end
    polyhedron = Polyhedron('A', A_poly, 'B', b_poly);
end


%% Comparing my elliptical approximation with maximal invariant set example from slides (Lec. 4 10-11)

A = [1 1;
    0 1];
B = [0; 1];
Q = eye(2);
R = 3;
[P,K,L] = idare(A,B,Q,R); K=-K;

F = [-1 0;
    0 1;
    zeros(2, 2)];
E = [zeros(2, 2);
    -1 0;
    0 1];
e = [ones(2, 1) * 3/4;
    ones(2, 1) * 1/2];

[Xf_H, Xf_h] = calculate_Xfa(A, B, K, P, F, E, e, [0; 0], [0; 0]);
my_inv_set = Polyhedron('A', Xf_H, 'B', Xf_h);


% computes a control invariant set for LTI system x^+ = A*x+B*u
system = LTISystem('A', A, 'B', B);
system.x.min = [-3/4, -3/4];
system.x.max = [3/4, 3/4];
system.u.min = [-1/2];
system.u.max = [1/2];
InvSet = system.invariantSet();

figure();
plot(my_inv_set, 'color', 'blue', 'alpha', 0.4); hold on;
plot(InvSet, 'color', 'red', 'alpha', 0.2);
axis equal
