%% Блок 1
f = @(x) sqrt(x) + cos(x);

a = 0;
b = 5;
h_x = 0.5;
h_xx = 0.001;
x = a : h_x : b;
xx = a : h_xx : b;

compareInterp(x, xx, f);

%%
g = @(x) x .^ 2 + cos(x); 
h = @(x) sin(20 * x);

max_derivative_g = 3;
max_derivative_h = 400;
a_priori_error_g = ones(size(xx, 2));
a_priori_error_h = ones(size(xx, 2));
j = 1;
w = (x(2) + x(1)) / 2;

a_priori_error_g = a_priori_error_g * max_derivative_g * abs((w - x(2)) * (w - x(1))) / 2;
a_priori_error_h = a_priori_error_h * max_derivative_h * abs((w - x(2)) * (w - x(1))) / 2;

figure(1)
hold on
title('График погрешности для функции g(x)');
plot(xx, a_priori_error_g, 'r');
a_posteriori_error_g = abs(g(xx) - interp1(x, g(x), xx, 'linear'));
plot(xx, a_posteriori_error_g, 'b');
hold off
legend('a priori', 'a posteriori');

figure(2)
hold on
title('График погрешности для функции h(x)');
plot(xx, a_priori_error_h, 'r');
a_posteriori_error_h = abs(h(xx) - interp1(x, h(x), xx, 'linear'));
plot(xx, a_posteriori_error_h, 'b');
hold off
legend('a priori', 'a posteriori');

%% Описания функций

function [] = compareInterp(x, xx, f)

    figure(1);
    hold on;
    title('NEAREST-интерполирование');
    plot(xx, f(xx), 'r');
    tic
    plot(xx, interp1(x, f(x), xx, 'nearest'), 'b');
    toc
    hold off; 
    legend('Красный - график исходной функции', 'Синий - NEAREST-интерполирование');
    
    figure(2);
    hold on;
    title('LINEAR-интерполирование');
    plot(xx, f(xx), 'r');
    tic
    plot(xx, interp1(x, f(x), xx, 'linear'), 'b');
    toc
    hold off; 
    legend('Красный - график исходной функции', 'Синий - LINEAR-интерполирование');

    figure(3);
    hold on;
    title('SPLINE-интерполирование');
    plot(xx, f(xx), 'r');
    tic
    plot(xx, interp1(x, f(x), xx, 'spline'), 'b');
    toc
    hold off; 
    legend('Красный - график исходной функции', 'Синий - SPLINE-интерполирование');
    
    figure(4);
    hold on;
    title('CUBIC-интерполирование');
    plot(xx, f(xx), 'r');
    tic
    plot(xx, interp1(x, f(x), xx, 'pchip'), 'b');
    toc
    hold off; 
    legend('Красный - график исходной функции', 'Синий - CUBIC-интерполирование');
end
