f = @(x) sqrt((x(1)- 0)^2 + (x(2) - 0)^2) - 1;

x0 = [0, 0];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];

rho = supportLebesgue(f, x0, A, b, Aeq, beq, lb, ub);       

% Центр ромба
xc_r = [2, 1];
% Половина первой диагонали
diag1 = 3;
% Половина второй диагонали
diag2 = 1.5;
% Ромб:
f_romb = @(x) (abs(x(1) - xc_r(1))/diag1) + (abs(x(2) - xc_r(2))/diag2) - 1;
% Функция rho для ромба:
rho_rhombus = @(l) [ max( diag1*abs(l(1)), diag2*abs(l(2)) ) + l(1)*xc_r(1) + l(2)*xc_r(2), xc_r(1) + (diag1*abs(l(1))>=diag2*abs(l(2)))*diag1*sign(l(1)), xc_r(2) + (diag1*abs(l(1))<diag2*abs(l(2)))*diag2*sign(l(2)) ];

% Центр квадрата
xc_s = [0, 0];
% Сторона квадрата
side_length = 1;
% Квадрат:
f_square = @(x) max(abs(x(1) - xc_s(1)), abs(x(2) - xc_s(2))) - side_length/2;
% Функция rho для квадрата:
rho_square = @(l) [ (side_length/2)*(abs(l(1)) + abs(l(2))) + l(1)*xc_s(1) + l(2)*xc_s(2), (l(1)>0)*side_length/2 + (l(1)<=0)*(-side_length/2)+xc_s(1), xc_s(2)+(l(2)>0)*side_length/2 + (l(2)<=0)*(-side_length/2) ]

drawPolar(rho_square, 80);

function [ supp_func_and_supp_vector ] = supportLebesgue(f, x0, A, b, Aeq, beq, lb, ub)
    scalar_product = @(x, l) x(1) * l(1) + x(2) * l(2);      
    function [c, ceq] = nlcon(x)
        c = f(x);
        ceq = [];
    end
    nonlcon = @nlcon;
    x_we_need = @(l) fmincon(@(x) -scalar_product(x, l), x0, A, b, Aeq, beq, lb, ub, nonlcon);
    supp_func_and_supp_vector = @(l) [ scalar_product(x_we_need(l), l), x_we_need(l)];
end



function [] = drawPolar(rho, N)
    array_for_polar = zeros(N, 2);
    for i = 0 : N - 1
        angle = 2*pi*i/N;
        next_angle = 2*pi*(i + 1)/N;
        
        vector_on_the_ring(1) = cos(angle);
        vector_on_the_ring(2) = sin(angle);
        next_vector(1) = cos(next_angle);
        next_vector(2) = sin(next_angle);
        
        value_and_point = rho(vector_on_the_ring);
        next_value_and_point = rho(next_vector);
        
        value = value_and_point(1);
        next_value = next_value_and_point(1);
        
        point = value_and_point(2:3);
        next_point = next_value_and_point(2:3);     
     
        hold on;
        line([point(1), next_point(1)], [point(2), next_point(2)], 'Color', 'blue');
        
        if value > 0
            array_for_polar(i + 1, 1) = vector_on_the_ring(1) / value; 
            array_for_polar(i + 1, 2) = vector_on_the_ring(2) / value;
            rho([array_for_polar(i+1,1), array_for_polar(i+1,2)]);
        elseif value < 0
            vector_on_the_ring = vector_on_the_ring .* (-1);
            value_and_point = rho(vector_on_the_ring);
            value = value_and_point(1);
            array_for_polar(i + 1, 1) = vector_on_the_ring(1) / value; 
            array_for_polar(i + 1, 2) = vector_on_the_ring(2) / value;
            rhooo = rho([array_for_polar(i+1,1), array_for_polar(i+1,2)]);
        end
    end
    [k, av] = convhull(array_for_polar);
    array = zeros(size(k, 1), 2);
    for i = 1:size(k, 1)
        array(i, 1) = array_for_polar(k(i), 1);
        array(i, 2) = array_for_polar(k(i), 2);
        
    end
    patch(array(:, 1), array(:, 2), 'red');
    
    axis([-2 8 -2 8]);
    hold off;
end

