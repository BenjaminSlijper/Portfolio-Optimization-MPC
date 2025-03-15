system = ltvss(@ltvStockModel, 1);

t = 0:11;
x0 = [0 0 0 0 10000]';
u = [ones(12,4)*100, zeros(12, 4)];

[y, ~, x] = lsim(system, u, t, x0);

plot(t, x)
