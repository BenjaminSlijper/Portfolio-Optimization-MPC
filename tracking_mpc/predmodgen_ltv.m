function [P,S] = predmodgen_ltv(ltvFunction,k,dim)
% Prediction matrices generation
%
% This function computes the prediction matrices P and S so that:
% y_N = P * x_0 + S * u_N
%
% P = A(k+N-1)A(k+N-2) ... A(k)
% S = [A(k+N-1)A(k+N-2)...B(k), A(k+N-2)...B(k+1), ... , B(k+N-1)]

% Prediction matrix from initial state

P = eye(size(A0));

for i=0:dim.N-1
    [A, ~] = ltvFunction(k+i);
    P = A * P;
end

% Prediction matrix from input

S = zeros(dim.nx, dim.N*dim.nu);
% Filling the array with [B(k), B(k+1), ..., B(k+N-1)]
for i=0:dim.N-1
    [~, B, ~] = ltvFunction(k+i);
    S(1:dim.nx, i*dim.nu+1: (i+1)*dim.nu) = B;
end

% Adding the A-matrices
for i=0:dim.N-2  % Loop over the blocks in the S-matrix
    cum_A = eye(dim.nx, dim.nx);
    start = dim.N - 1;
    ending = i + 1;
    
    % Building the A(k+N-1) * A(k+N-2) * ... matrix
    for j=start:-1:ending
        [A, ~] = ltvFunction(k+j);
        disp(k+j)
        cum_A = A * cum_A;
    end
    S(1:dim.nx, i*dim.nu+1: (i+1)*dim.nu) = cum_A * S(1:dim.nx, i*dim.nu+1: (i+1)*dim.nu);
end