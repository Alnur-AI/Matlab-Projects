rad = 5; % радиус окружности
t_end = 154.78;
tspan = [0, t_end]; % промежуток времени
v = [11, -5]; % начальная скорость
alpha = 2.54;
start_point = [3, -0.2987];



options = odeset('Events', @reached_the_circle_Event);
[t, x] = ode45( @(t, x) [v(1); v(2)], tspan, start_point, options );
draw_dot(t, x, rad);
while t(end) ~= t_end
    t(end);
    v_before_change = v
    if x(length(t), 1) ~= 0 & x(length(t), 2) ~= 0
        k = x(length(t), 1)
        l = x(length(t), 2)
        c = -v(1)*k - v(2)*l;
        z = c / k;
        a = (k^2+l^2)/(k^2);
        b = (-2)*z*l/k;
        d = z^2 - v(1)^2 - v(2)^2;
        D = b^2 - 4 * a * d;
        if D < 0 
            D = round(D);
        end
        if abs((-b + sqrt(D))/(2*a) + v(2)) < 0.01
            v2 = (-b - sqrt(D))/(2*a);
            v1 = z - (l/k)*v2;
        else
            v2 = (-b + sqrt(D))/(2*a)
            v1 = z - (l/k)*v2
        end
    elseif x(length(t), 1) == 0
        
        v1 = -v(1);
        v2 = v(2);
    elseif x(length(t), 2) == 0
        
        v1 = v(1);
        v2 = -v(2);
    end
    v_after_change = [v1, v2]
    v = [v1, v2] ./ alpha;
    tspan = [t(end), t_end];
    eps = 0.05;
    [t, x] = ode45( @(t, x) [v(1); v(2)], tspan, [x(length(t), 1); x(length(t), 2)] + [v(1); v(2)] .* eps, options );
    draw_dot(t, x, rad);
end

function [value, isterminal, direction] = reached_the_circle_Event(t, x)
    rad = 5;
    value = x(1)^2 + x(2)^2 - rad^2;
    isterminal = 1;   
    direction = 0;   
end

function [] = draw_dot(t, x, rad)
    grid = -rad : 0.01 : rad;
    for i = 1:length(t)
        plot(grid, sqrt(ones(1, length(grid)) .* rad^2 - grid.^2), 'b');
        hold on;
        plot(grid, -sqrt(ones(1, length(grid)) .* rad^2 - grid.^2), 'b');
        axis([ - rad - 1, rad + 1, - rad - 1, rad + 1]);
        plot(x(i, 1), x(i, 2), '-o', 'MarkerIndices', 1, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'green', 'MarkerSize', 4);
        plot(x(1:i, 1), x(1:i, 2)), '-k';
        hold off;
        pause(0.01);
    end
    
end
