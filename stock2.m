%% Пункты 1 и 2
clear;
clc;

% Объем выборки случайной величины X (здесь и далее X - канторова случайная величина)
n = 100;
% Объем выборки случайной величины Y = 1 - X
m = 90; 

alpha = 0.1; % Уровень значимости
times = 1000; % Количество испытаний

counter_cant = 0; % Счетчик принятия гипотезы о том, что случайная величина X является канторовой
counter_symmetry = 0; % Счетчик принятия гипотезы о том, что случайные величины X и Y=1-X одинаково распределены
counter_self_similarity = 0; % Счетчик принятия гипотезы о том, случайная величина X/3 распределена так же, как и с. в. X при условии X<=1/3


for i = 1:times
    % Генерируем выборку и вектор значений функции распределения случайной величины X
    [x, F] = cantrnd(n);
    % Вычисляем статистику Колмогорова 
    empirF = sum((ones(n, n) .* x' < x), 2)./n;
    Dn = max(abs(empirF - F));
    % Принимаем либо отвергаем гипотезу
    p = 1 - F_K(sqrt(n) * Dn);
    counter_cant = counter_cant + (p > alpha);
    
    % Генерируем выборку случайной величины Y = 1 - X
    y = cantrnd(m);
    y = ones(m, 1) - y;
    % Вычисляем статистику Смирнова
    empF_minus_empG = (sum((ones(n + m, n).* x' < cat(1,x,y)), 2) ./ n) - (sum((ones(n + m, m).* y' < cat(1,x,y)), 2) ./ m);
    Dnm = max(abs(empF_minus_empG));
    % Принимаем либо отвергаем гипотезу
    p = 1 - F_K(sqrt(n*m/(n + m)) * Dnm);
    counter_symmetry = counter_symmetry + (p > alpha);
    
    % Генерируем выборки случайных величин X при условии X <= 1/3 и X/3
    m = sum(x <= 1/3); % объем выборки случайной величины X при условии X <= 1/3
    y = sort(x.*(x <= 1/3));
    y = y(n - m + 1 : end);
    x = x / 3;
    % Вычисляем статистику Смирнова
    empF_minus_empG = (sum((ones(n + m, n).* x' < cat(1,x,y)), 2) ./ n) - (sum((ones(n + m, m).* y' < cat(1,x,y)), 2) ./ m);
    Dnm = max(abs(empF_minus_empG));
    % Принимаем либо отвергаем гипотезу
    p = 1 - F_K(sqrt(n*m/(n + m)) * Dnm);
    counter_self_similarity = counter_self_similarity + (p > alpha);
end

str=sprintf('Частота принятия гипотезы о том, что выборка сгенерирована при помощи сингулярного канторова распределения, равна:');
disp(str);
counter_cant/times   

str=sprintf('Частота принятия гипотезы о том, что случайные величины X и Y = 1 - X одинаково распределены (здесь X - канторова с. в.):');
disp(str);
counter_symmetry/times

str=sprintf('Частота принятия гипотезы о том, что канторова случайная величина X обладает свойством самоподобия относительно деления на 3:');
disp(str);
counter_self_similarity/times

%% Картинки к пунктам 1, 2

% Иллюстрируем Критерий Колмогорова
num = 10000;
[x, Fx] = cantrnd(num);
array_x = sortrows(cat(2, x, Fx));

arg_x = array_x(:,1);
expectedF_x = array_x(:,2);
empiricalF_x = sum(ones(num, num) .* arg_x' < arg_x, 2) / num;

figure;
hold on;
plot(arg_x, empiricalF_x)
plot(arg_x, expectedF_x)
legend('Эмпирическая функция распределения', 'Теоретическая функция распределения');
hold off; 

% Иллюстрируем свойство симметричности относительно 0.5
y = cantrnd(num);
y = ones(num, 1) - y;
y = sort(y);

empiricalF_y = sum(ones(num, num) .* y' < y, 2) / num;

figure;
hold on;
plot(arg_x, empiricalF_x)
plot(y, empiricalF_y)
legend('Эмпирическая функция распределения с. в. X', 'Эмпирическая функция распределения с. в. 1-X');
hold off;

% Иллюстрируем свойство самоподобия
x = cantrnd(num);
m = sum(x <= 1/3);
y = sort(x.*(x <= 1/3));
y = y(num - m + 1 : end);
x = sort(x) / 3;

empiricalF_y = sum(ones(m, m) .* y' < y, 2) / m;
empiricalF_x = sum(ones(num, num) .* x' < x, 2) / num;

figure;
hold on;
plot(x, empiricalF_x)
plot(y, empiricalF_y)
legend('Эмпирическая функция распределения с. в. X/3', 'Эмпирическая функция распределения X при условии X\leq1/3');
hold off;
%% Пункт 3

N =1000;

sample_mean = zeros(N, 1);
sample_dispersion = zeros(N, 1);
for i = 1:N
    x = cantrnd(i);
    sample_mean(i) = (1/i) * sum(x);
    sample_dispersion(i) = (1/i) * sum( (x - ones(i,1).*sample_mean(i)).^2 );
end

figure(1);
hold on;
loglog(1:N, sample_mean);
plot(1:N, ones(N,1)*0.5);
legend('Выборочное среднее', 'Математическое ожидание');
hold off;

figure(2);
hold on;
plot(1:N, sample_dispersion);
plot(1:N, ones(N,1)*0.125);
legend('Выборочная дисперсия', 'Дисперсия');
hold off;

%% Функции
% Функция распределения Колмогорова
function result = F_K(x, n)
    if nargin < 2
        n = 1000;
    end
    result = 1;
    for k = 1:n
        result = result + 2*((-1)^k)*exp(-2*(k^2)*(x^2));
    end
end

% Функция генерации канторовых случайных величин
function [x, F] = cantrnd(n, eps)
    if nargin < 1
        n = 1; % объем выборки
    end
    if nargin < 2
        eps = 1e-10; % точность расчета
    end
    m = round(-log(eps)/log(3));
    bern = bernrand(0.5, n, m);
    deg = -(1:m)';
    x = 2 * bern * 3.^deg;
    F = bern * 2.^deg;
end

% Функция генерации случайных величин Бернулли
function x = bernrand(p, varargin)
    if (p>1) || (p<0)
        error('p = %4.2f is out of bounds [0, 1]',p)
    end
    x = rand(varargin{:}) < p;
end