function [T,S,p] = predmodgen_ltv(ltvFunction,k,dim)
    T = zeros(dim.nx * (dim.N + 1), dim.nx);
    S = zeros(dim.nx * (dim.N + 1), dim.nu * dim.N);

    power_matrices = {eye(dim.nx)};
    for i=1:dim.N
        [A, ~] = ltvFunction(k+i-1);
        power_matrices{end+1} = power_matrices{i} * A;
    end

    for i=0:dim.N
        T(i * dim.nx + 1: (i+1) * dim.nx, :) = power_matrices{i+1};
        for j=0:dim.N
            if i > j
                [~, B, ~] = ltvFunction(k+j);
                S(i * dim.nx + 1: (i+1) * dim.nx, j * dim.nu + 1: (j+1) * dim.nu) = power_matrices{i - j} * B;
            end
        end
    end
    p = power_matrices;
end