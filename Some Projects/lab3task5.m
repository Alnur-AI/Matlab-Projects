a = 4; % альфа
b = 1; % бета
g = 10; % гамма
d = 0.23; % дельта
x0 = 5;
y0 = 1;

tspan = [0, 10]; % интервал интегрирования 

[t, x_y] = ode45(@(t, x_y) [ a*x_y(1) - g*x_y(1)*x_y(2); -b*x_y(2) + d*x_y(1)*x_y(2) ], tspan, [x0; y0]);
hold on;
plot(t, x_y(:, 1), 'b');
plot(t, x_y(:, 2), 'r');
legend('x(t)', 'y(t)')
hold off;
figure(2);
plot3(x_y(:, 1), x_y(:, 2), t);



% аналитическое решение для параметров a=b=g=1, d=0, x0=y0=1
figure(3)
grid1 = 0:0.01:10;
hold on; 
title('Аналитическое решение');
plot(grid1, exp(grid1 + exp(-grid1) - 1), 'r');
plot(grid1, exp(-grid1), 'b');
legend('x(t)', 'y(t)')
hold off;

figure(4)
title('');
plot3(exp(grid1 + exp(-grid1) - 1), exp(-grid1), grid1);