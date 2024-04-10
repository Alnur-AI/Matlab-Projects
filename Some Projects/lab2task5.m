f = @(x) exp(x ./ 2) - x .^ 3 - x;

a = -1;
b = 1;
n = 300;

fourierApprox (f, a, b, n, 'lezhandr');

function [] = fourierApprox (f, a, b, n, meth)
    
    if  ~strcmp(meth, 'trig') & ~strcmp(meth, 'cheb') & ~strcmp(meth, 'lezhandr')
        error('Ошибка! Некорректно указана система функций')
    end

    
    mov(1 : n) = struct('cdata', [], 'colormap', []); % Выделим память под кадры анимации
    h = 0.001;
    x = a : h : b;
    sum = zeros(1, size(x, 2)); % Частичная сумма ряда Фурье
    
    
    
    for i = 1:n
        plot(x, f(x), 'r');
        axis([-1 1 -0.5 3]);
        hold on;
        if strcmp(meth, 'trig')
            left = -pi;
            right = pi;
            func_i = getFunc_trig(i);
            weight_func = @(x) ones(1, size(x, 2));
        elseif strcmp(meth, 'cheb')
            left = -1;
            right = 1;
            func_i = getFunc_cheb(i);
            weight_func = @(x) ones(1, size(x, 2)) ./ sqrt(1 - x .^ 2);
        elseif strcmp(meth, 'lezhandr')
            left = -1;
            right = 1;
            func_i = getFunc_lezhandr(i);
            func_i = @(x) func_i(x) ./ sqrt(2 / (2 * i - 1));
            weight_func = @(x) ones(1, size(x, 2));
        end
        
        func_find_lyambda= @(x) (func_i(x) .^ 2) .* weight_func(x);
        func_find_c = @(x) f(x) .* func_i(x) .* weight_func(x);
        lyambda_i = integral(func_find_lyambda, left, right);
        c_i = integral(func_find_c, left, right) / lyambda_i;
        sum = sum + func_i(x) .* c_i;
        
        plot(x, sum, 'b');
        hold off;
        title(['Частичная сумма ряда Фурье по', newline, 'системе ', meth, ' с номером ', num2str(i)]);
        mov(i) = getframe();
        pause(0.1);
    end
end

function [ f ] =  getFunc_trig(i)
    if i == 1
        f = @(x) ones(1, size(x, 2)) / sqrt(2 * pi);
    elseif mod(i, 2) == 1
        f = @(x) sin(floor(i / 2) * x) / sqrt(pi);
    else
        f = @(x) cos((i / 2) * x) / sqrt(pi);
    end
end

function [ f ] = getFunc_cheb(i)
    if i == 1
        f = @(x) cos((i - 1) * acos(x)) ./ sqrt(pi);
    else
        f = @(x) cos((i - 1) * acos(x)) ./ sqrt(pi / 2);
    end
end


function [ f ] =  getFunc_lezhandr(i)
    if i == 1
        f = @(x) ones(1, size(x,2));
    elseif i == 2
        f = @(x) x;
    else
        lezh1 = getFunc_lezhandr(i - 1);
        lezh2 = getFunc_lezhandr(i - 2);
        f = @(x)  x .* (2 * i - 3) ./ (i - 1) .* lezh1(x) - lezh2(x) .* (i - 2) ./ (i - 1); 
    end
end

