function [A,B,C,D,E,dx0,x0,u0,y0,Delays] = ltvStockModel(k)

    % Declare persistent variables
    persistent apple boeing nvidia oneok;
    
    % If data is not yet loaded, read from CSV files
    if isempty(apple) || isempty(boeing) || isempty(nvidia) || isempty(oneok)
        % If paths do not work, run init.m to add all files to the path
        apple = readmatrix("data_processed\apple_return.csv");
        boeing = readmatrix("data_processed\boeing_return.csv");
        nvidia = readmatrix("data_processed\nvidia_return.csv");
        oneok = readmatrix("data_processed\oneok_return.csv");
    end
    
    interest = 0.000274; 
    alpha = 0.0001; 
    beta = 0.0001;
    
    apple_return = apple(:, 2);
    boeing_return = boeing(:, 2);
    nvidia_return = nvidia(:, 2);
    oneok_return = oneok(:, 2);


    
    % Plot the Apple return data (only the second column as it's assumed to be the return values)
    % Index is k+1 instead of k due to MATLAB being 1-indexed
    A = diag([1 + apple_return(k+1); ...
              1 + boeing_return(k+1); ...
              1 + nvidia_return(k+1); ...
              1 + oneok_return(k+1); ...
              1 + interest]);
    
    return_matrix = diag([1 + apple_return(k+1); ...
              1 + boeing_return(k+1); ...
              1 + nvidia_return(k+1); ...
              1 + oneok_return(k+1)]);
    
    final_row_B_1 = ones(1, length(return_matrix)) * (1 + interest) * (-1 - alpha);
    final_row_B_2 = ones(1, length(return_matrix)) * (1 + interest) * (1 - beta);
    
    B = [return_matrix -return_matrix;
        final_row_B_1 final_row_B_2];
    
    C = ones(1, length(A));

    D = zeros(1, length(B(1, :)));
    
    E = [];
    dx0 = [];
    x0 = [];
    u0 = [];
    y0 = [];
    Delays = [];
end