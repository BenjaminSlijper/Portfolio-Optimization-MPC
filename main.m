dT = 1;

system = ltvss(@ltvStockModel, "Ts", dT, "TimeUnit", 'days');

t = 0:dT:11;
x0 = [0 0 0 0 10000]';
u = [ones(5,4)*1000, zeros(5, 4) ; zeros(7, 8)];

figure()
lsim(system, u, t, x0)
hold on
[y,t,x,p] = lsim(system, u, t, x0);
plot(t, x)
hold off
