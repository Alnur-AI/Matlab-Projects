% Настройка параметров системы
% Пример 1
A = [ 1, 2; -2, -1 ];
B = [ 5, 4; -8, 1 ];
f = [4; 0];
t0 = 0;
p = [2, 1];
alpha = 3;
beta = 0.8;
gamma = 5;
k = 3;
x0 = [-1; 1];
x1 = [2.2; 0.8];
% Пример 2
%A = [ 1, 1; 3, 1 ];
%B = [ 0, 1; 0, 11 ];
%f = [5; 5];
%t0 = 0;
%p = [4, 1];
%alpha = 1;
%beta = 2;
%gamma = 3;
%k = 2;
%x0 = [2; 5];
%x1 = [6.6816, 8.1905];

% Рисуем начальное множество Х0
step = 0.01;
[X, Y] = meshgrid(  x0(1) - abs(k) : step : x0(1) + abs(k),  x0(2) - abs(k) : step : x0(2) + abs(k)  );
Z = zeros(size(Y, 1), size(Y, 2));
for i = 1:size(Z, 1)
    for j = 1:size(Z, 2)
        if (abs(X(i,j) - x0(1)) <= k/2) && (abs(Y(i,j) - x0(2)) <= k/2) 
            Z(i, j) = 1;
        end
    end
end
hold on;
contourf(X, Y, Z, [0.5, 0.5]);

% Настройка параметров численного метода 
num = 500;
h = 2 * pi / num;
t1 = 0.5; % выберем параметр t1 достаточно большим и будет интегрировать дифф. ур. d(psi)/dt=-A'psi по интервалу [t0, t1]
T = t1 - t0;
iter = 1;

% Рисуем основную картинку
for i = 1 : (num - 1)
    [T, iter] = Draw_main_picture(i, h, A, B, f, t0, t1, k, x0, p, alpha, beta, gamma, x1, T, iter, 0);
end

% Выводим на экран оптимальное время
T

% Отмечаем зеленым маркером целевую точку
plot(x1(1), x1(2), '-ob', 'MarkerIndices', 1, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'green', 'MarkerSize', 5.5);

% Выделяем отдельным цветом оптимальную траекторию 
[T, iter, psi_t, optpsi, time, optx, optu] = Draw_main_picture(iter, h, A, B, f, t0, t1, k, x0, p, alpha, beta, gamma, x1, T, iter, 1);

% Зависимость первой компоненты оптимальной траектории от времени
figure;
plot(time, optx(:, 1));
%title('Зависимость первой компоненты оптимальной траектории от времени');

% Зависимость второй компоненты оптимальной траектории от времени
figure;
plot(time, optx(:, 2));
%title('Зависимость второй компоненты оптимальной траектории от времени');

% Зависимость первой компоненты оптимального управления от времени
figure;
plot(psi_t, optu(:, 1));
%title('Зависимость первой компоненты оптимального управления от времени');

% Зависимость второй компоненты оптимального управления от времени
figure;
plot(psi_t, optu(:, 2));
%title('Зависимость второй компоненты оптимального управления от времени');

% Зависимость первой компоненты сопряженных переменных от времени 
figure;
plot(psi_t, optpsi(:, 1));
%title('Зависимость первой компоненты сопряженных переменных от времени');

% Зависимость второй компоненты сопряженных переменных от времени
figure;
plot(psi_t, optpsi(:, 2));      
%title('Зависимость второй компоненты сопряженных переменных от времени');

hold off;

%err = (scal_product(-optpsi(end, :), optx(end, :)) - scal_product(-optpsi(end, :), x1)) / sqrt(optpsi(end,1)^2 + optpsi(end,2)^2)

function answ = scal_product(l1, l2)
    answ = l1(1)*l2(1)+l1(2)*l2(2);
end

function [T, iter, psi_t, optpsi, time, optx, optu] = Draw_main_picture (i, h, A, B, f, t0, t1, k, x0, p, alpha, beta, gamma, x1, curr_T, curr_iter, optimal)
    psi_t0 = [cos(h * i), sin(h * i)];
    [t, psi] = ode45(@(t, psi) fun_ATpsi(A, t, psi), [t0, t1], psi_t0);
    psi_t = t;
    optpsi = psi;
    % Из второго условия в ПМП находим x(t0):
    x_t0 = [ x0(1) + (k/2) * sign(psi_t0(1)),  x0(2) + (k/2) * sign(psi_t0(2)) ];
    % Находим "кандидата" на оптимальное управление u*
    u_star = zeros(size(psi, 1), 2);
    for j = 1:size(psi, 1)
        BTpsi_j = B' * psi(j, :)';
        if BTpsi_j(1) < 0 
            u_star(j, :) = [p(1) + (BTpsi_j(1) / alpha) * sqrt(alpha * gamma / (BTpsi_j(1)^2 + BTpsi_j(2)^2)), p(2) + (BTpsi_j(2) / alpha) * sqrt(alpha * gamma / (BTpsi_j(1)^2 + BTpsi_j(2)^2))];
        else
            u_star(j, :) = [p(1) + (BTpsi_j(1) / beta) * sqrt(alpha * beta * gamma / (alpha*BTpsi_j(1)^2 + beta*BTpsi_j(2)^2)), p(2) + (BTpsi_j(2) / alpha) * sqrt(alpha * beta * gamma / (alpha*BTpsi_j(1)^2 + beta*BTpsi_j(2)^2))];
        end
    end
    
    % Решаем основной диффур
    options = odeset('Events', @(t, x) reached_X1_event(t, x, x1));
    [t_in_main_diff, x] = ode45(@(t_in_main_diff, x) fun_main_diff(A, B, f, t, u_star, t_in_main_diff, x), [t0, t1], x_t0, options);
    for l = 1:length(t_in_main_diff)
        if sqrt((x(1)-x1(1))^2 + (x(2)-x1(2))^2) < 0.05
            x = x(1:l, :);
            break;
        end
    end
    
    % Рисуем i-ую траекторию 
    if optimal == 0
        %plot(x(:, 1), x(:, 2), '-ob', 'MarkerIndices', 1:1:size(x,1), 'MarkerEdgeColor', 'yellow', 'MarkerFaceColor', 'yellow', 'MarkerSize', 2);
        plot(x(:, 1), x(:, 2), 'b');
    else
        plot(x(:, 1), x(:, 2), 'r');
    end
    if (t_in_main_diff(end) - t0 < curr_T)
        T = t_in_main_diff(end) - t0;
        iter = i;
    else
        T = curr_T;
        iter = curr_iter;
    end
    
    time = t_in_main_diff;
    optx = x;
    optu = u_star;
    
end

function result = fun_ATpsi(A, t, psi)
    result = [-(A(1, 1) * psi(1) + A(2, 1) * psi(2)),  -(A(1, 2) * psi(1) + A(2, 2) * psi(2))]';
end

function [value, isterminal, direction] = reached_X1_event(t, x, x1)
    r = sqrt((x1(1) - x(1))^2 + (x1(2) - x(2))^2);
    if r < 0.05
        value = 0;
    else
        value = 1;
    end
    isterminal = 1;   
    direction = 0;   
end

function dotx = fun_main_diff(A, B, f, t, u, t_in_md, x)
    dotx = zeros(2, 1);
    dotx(1) = x(1)*A(1,1) + x(2)*A(1,2) + B(1,1)*interp1(t, u(:, 1), t_in_md) + B(1, 2)*interp1(t, u(:, 2), t_in_md) + f(1); 
    dotx(2) = x(1)*A(2,1) + x(2)*A(2,2) + B(2,1)*interp1(t, u(:, 1), t_in_md) + B(2, 2)*interp1(t, u(:, 2), t_in_md) + f(2);
end




