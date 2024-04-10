figure_for_plotting = gcf;

plotFT(figure_for_plotting, @func2, @ftfunc2, 0.001, [-100 100], [-10 10]);

function [ structure ] = plotFT (hFigure, fHandle, fFTHandle, step, inpLimVec, outLimVec)
    SPlotInfo = get(hFigure, 'UserData');
    if ~isempty(SPlotInfo)
        delete(SPlotInfo.first_axises);
        delete(SPlotInfo.second_axises);
    else
        SPlotInfo = struct('first_axises', [], 'second_axises', []);
    end
    
    a = inpLimVec(1);
    b = inpLimVec(2);
    T = b - a;  
    N = T / step;
    N = round(N);
    step = T / N;
    grid_ab = a : step : b;
    func_on_grid_ab = fHandle(grid_ab);
    
    n = 0;
    if a > 0 && b > 0
        while n * T < a
            n = n + 1;
        end
    elseif a < 0 && b < 0
        while n * T > b 
            n = n - 1;
        end
    end
    vodorazdel = n * T;
    
    idx = 1;
    while grid_ab(idx + 1) <= vodorazdel
        idx = idx + 1;
    end
    %grid_0T = 0 : step : T;
    %func_on_grid_0T = zeros(1, N + 1);
    func_on_grid_0T(1 : (N + 1) - idx) = func_on_grid_ab(idx + 1 : N + 1);
    func_on_grid_0T((N + 1) - idx + 1 : N + 1) = func_on_grid_ab(1 : idx);
    
    length(func_on_grid_0T);
    length(func_on_grid_ab);
    
    fourier = step * fft(func_on_grid_0T);
            
    new_step = (2 * pi) / T;
    new_grid = 0 : new_step : new_step * N;
    
    c = outLimVec(1);
    d = outLimVec(2);
    
    new_T = new_step * N
    
    left = -new_T;
    right = new_T;
    
    counter = 2;
    while c < left || d > right
        left = left - new_T;
        right = right + new_T;
        counter = counter + 2;
    end
    big_grid = left : new_step : right;
    fourier_on_big_gr = repmat(fourier(2:end), 1, counter);
    length(big_grid(2:end));
    
    subplot(2, 1, 1);
    hold on;
    plot(big_grid(2:end), real(fourier_on_big_gr), 'b');
    %plot(grid_0T, func_on_grid_0T, 'r');
    legend('Вещественная часть аппроксимации F(\lambda) через БПФ');
    if ~isempty(fFTHandle)
        grid_cd = outLimVec(1) : 0.01 : outLimVec(2);
        plot(grid_cd, real(fFTHandle(grid_cd)) , 'r', 'DisplayName', 'Вещественная часть аналитически посчитанного F(\lambda)');
    end
    axis([outLimVec(1) outLimVec(2) min(real(fourier_on_big_gr)) max(real(fourier_on_big_gr))]);
    ylabel('Re F(\lambda)');
    xlabel('\lambda');
    SPlotInfo.first_axises = gca;
    hold off;   
    
    subplot(2, 1, 2);
    hold on;
    plot(big_grid(2:end), imag(fourier_on_big_gr), 'b');
    %plot(grid_0T, func_on_grid_0T, 'r');
    legend('Мнимая часть аппроксимации F(\lambda) через БПФ');
    if ~isempty(fFTHandle)
        grid_cd = outLimVec(1) : 0.01 : outLimVec(2);
        plot(grid_cd, imag(fFTHandle(grid_cd)) , 'r', 'DisplayName', 'Мнимая часть аналитически посчитанного F(\lambda)');
    end
    axis([outLimVec(1) outLimVec(2) min(imag(fourier_on_big_gr)) max(imag(fourier_on_big_gr))]);
    ylabel('Im F(\lambda)');
    xlabel('\lambda');
    SPlotInfo.second_axises = gca;
    hold off;
    
    structure = struct('nPoints', N, 'step', step, 'inpLimVec', inpLimVec, 'outLimVec', outLimVec);
    set(hFigure, 'UserData', SPlotInfo);
    
end

function func_value = func1(t)
    func_value = t .* exp(-t.^2);
end

function func_value = func2(t)
    for i = 1:length(t)
        if t(i) ~= 0
            func_value(i) = (cos(t(i)) - exp(-abs(t(i)))) ./ t(i);
        else
            func_value(i) = 0;
        end
    end
end

function func_value = func3(t)
    func_value = exp(-2 * abs(t)) ./ (1 + (cos(t)).^2);
end

function func_value = func4(t)
    func_value = 2 * ones(1, length(t)) ./ (3 + 4 * t.^4);
end

function func_value = ftfunc1(l)
    Re_part = zeros(1, length(l));
    Im_part = -sqrt(pi)/2 * l .* exp((-l.^2) / 4);
    func_value = Re_part + j * Im_part;
end

function func_value = ftfunc2(l)
    Re_part = zeros(1, length(l));
    Im_part = (pi*sign(-l - 1) + pi*sign(1 - l) - 4*atan(-l)) / 2;
    func_value = Re_part + j * Im_part;
end

function func_value = ftest(t)
    %func_value = ones(1, length(t)) ./ t;
    %func_value = ones(1, length(t));
    %for i = 1:length(t)
    %    if (t(i) <= 9) & (t(i) >= 2.5)
    %        func_value(i) = 1;
    %    else
    %         func_value(i) = 0;
    %    end
    %end
end