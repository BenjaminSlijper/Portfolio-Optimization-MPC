function [A,B,C,D,E,dx0,x0,u0,y0,Delays] = ltvStockModel(k)

    apple = readmatrix("data_processed\apple_return.csv");
    boeing = readmatrix("data_processed\boeing_return.csv");
    nvidia = readmatrix("data_processed\nvidia_return.csv");
    oneok = readmatrix(	"data_processed\oneok_return.csv");
    
    bank_return = 0.1; % store somewhere else
    alpha = 0.0001;  % store somewhere else
    beta = 0.0001;  % store somewhere else
    
    apple_return = apple(:, 2);
    boeing_return = boeing(:, 2);
    nvidia_return = nvidia(:, 2);
    oneok_return = oneok(:, 2);
    
    A = diag([1 + apple_return(k+1); ...
              1 + boeing_return(k+1); ...
              1 + nvidia_return(k+1); ...
              1 + oneok_return(k+1); ...
              1 + bank_return]);
    
    return_matrix = diag([1 + apple_return(k+1); ...
              1 + boeing_return(k+1); ...
              1 + nvidia_return(k+1); ...
              1 + oneok_return(k+1)]);
    
    final_row_B_1 = ones(1, length(return_matrix)) * (1 + bank_return) * (-1 - alpha);
    final_row_B_2 = ones(1, length(return_matrix)) * (1 + bank_return) * (1 - beta);
    
    B = [return_matrix -return_matrix;
        final_row_B_1 final_row_B_2];
    
    C = ones(1, length(A));
    % C(end) = -1;

    D = zeros(1, length(B(1, :)));
    E = [];
    dx0 = [];
    x0 = [];
    u0 = [];
    y0 = [];
    Delays = [];
end