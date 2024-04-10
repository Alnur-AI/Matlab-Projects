A = [3, 1; 1, 3]; % узел
B = [1, 0; 0, 1]; % дикритический узел
C = [-1, 0; 0, 1]; % седло
D = [0, 4; -1, -1]; % фокус
E = [0, 1; -0.012, 0]; % центр

hold on;

tspan = [-25, 25];

h = 0.1;
eps = 11;

for i = 1:2
    for j = 1:2
        m(i, j) = E(i, j);
    end
end

if m == A | m == B
    eps = 0.5;
end

for i = 0 : h : 2 * pi
    [t, y] = ode45(@(t, y) [ m(1, 1) * y(1) + m(1, 2) * y(2); m(2, 1) * y(1) + m(2, 2) * y(2) ], tspan, [sin(i); cos(i)] .* eps);
    plot(y(:, 1), y(:, 2), 'b');
    quiver(y(10, 1), y(10, 2), y(15, 1) - y(10, 1), y(15, 2) - y(10, 2), 0.1, 'k', 'MaxHeadSize', 10);
    grid on;
end

axis([-10 10 -10 10]);

hold off;











