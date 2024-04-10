% Поточечная и среднеквадратичная сходимость, но не равномерная 
f1 = @(x) x .* (x == 1) + 0 * (x < 1);
f1n = @(x, n) x .^ n;

% Равномерная сходимость (и поточечная, и среднеквадратичная)
f2 = @(x) 0;
f2n = @(x, n) sin(x) .^ n;

% Сходится поточечно, но не в среднем и не равномерно
f3 = @(x) 0;
f3n = @(x, n) n * (0 < x & x < 1/n) ;

a = 0;
b = 1;
n = 10;

convergenceFunc(f2n, f2, a, b, n, 'uniform'); % pointwise - поточечная, uniform - равномерная, RMS - среднеквадратичная

function [] = convergenceFunc (fn, f, a, b, n, convType)

    if  ~strcmp(convType, 'pointwise') & ~strcmp(convType, 'uniform') & ~strcmp(convType, 'RMS')
        error('Ошибка! Некорректно указан тип сходимости')
    end
    
    mov(1 : n) = struct('cdata', [], 'colormap', []); % Выделим память под кадры анимации
    h = 0.01;
    x = a : h : b;
    
    for i = 1:n    
        
        plot(x, f(x), 'r');
        hold on;
        
        if strcmp(convType, 'uniform')
            str = num2str(  max(abs(f(x) - fn(x, i))), '%.4f'  );
            title( ['Значение метрики разности', newline 'для равномерной сходимости: ',  str] );
        end
        
        if strcmp(convType, 'RMS')
            h = f(x) - fn(x, i);
            h = h .^ 2;
            str = num2str(  sqrt(trapz(x, h)), '%.4f'  );
            title( ['Значение метрики разности для', newline, 'среднеквадратичной сходимости: ',  str] );
        end
                
        plot(x, fn(x, i), 'b');
        hold off;
        
        mov(i) = getframe();
        pause(0.4);
    end
    %movie(mov);
end