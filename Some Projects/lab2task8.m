f = @(x) sqrt(x(1)^2 + x(2)^2) - 2;

% Центр ромба
xc_r = [2, 1];
% Половина первой диагонали
diag1 = 3;
% Половина второй диагонали
diag2 = 1.5;
% Ромб:
f_romb = @(x) (abs(x(1) - xc_r(1))/diag1) + (abs(x(2) - xc_r(2))/diag2) - 1;
% Функция rho для ромба:
supp_function_rhombus = @(l) [ max( diag1*abs(l(1)), diag2*abs(l(2)) ) + l(1)*xc_r(1) + l(2)*xc_r(2), xc_r(1) + (diag1*abs(l(1))>=diag2*abs(l(2)))*diag1*sign(l(1)), xc_r(2) + (diag1*abs(l(1))<diag2*abs(l(2)))*diag2*sign(l(2)) ];

% Центр квадрата
xc_s = [0, 0];
% Сторона квадрата
side_length = 2;
% Квадрат:
f_square = @(x) max(abs(x(1) - xc_s(1)), abs(x(2) - xc_s(2))) - side_length/2;
% Функция rho для квадрата:
supp_function_square = @(l) [ (side_length/2)*(abs(l(1)) + abs(l(2))) + l(1)*xc_s(1) + l(2)*xc_s(2), (l(1)>0)*side_length/2 + (l(1)<=0)*(-side_length/2)+xc_s(1), xc_s(2)+(l(2)>0)*side_length/2 + (l(2)<=0)*(-side_length/2) ]

% Центр эллипса
xc_e = [1, 1];
% Полуоси эллипса
a = 2;
b = 6.56;
% Эллипс:
f_ellipse = @(x) (x(1) - xc_e(1))^2 / a^2 + (x(2) - xc_e(2))^2 / b^2 - 1;
% Функция rho для эллипса
supp_function_ellipse = @(l) [ l(1)*xc_e(1) + l(2)*xc_e(2) + sqrt(  (l(1)*a)^2 + (l(2)*b)^2  ), xc_e(1)+l(1)*(a^2)/(sqrt((l(1)^2)*(a^2) + (l(2)^2)*(b^2))), xc_e(2)+l(2)*(b^2)/(sqrt((l(1)^2)*(a^2) + (l(2)^2)*(b^2)))]


drawSet(supp_function_ellipse, 26)

supp_function_square([1, 0])

function [] = drawSet(rho, N)
    for i = 0 : N - 1
        angle = 2*pi*i/N;
        next_angle = 2*pi*(i + 1)/N;
        vector_on_the_ring(1) = cos(angle);
        vector_on_the_ring(2) = sin(angle);
        next_vector(1) = cos(next_angle);
        next_vector(2) = sin(next_angle);
        
        segment_left = -8;
        segment_right = 8;
        h = 0.01;
        setka = segment_left : h : segment_right;
        
        value_and_point = rho(vector_on_the_ring);
        value = value_and_point(1);
        point = value_and_point(2:3);
        
        next_value_and_point = rho(next_vector);
        next_value = next_value_and_point(1);
        next_point = next_value_and_point(2:3);
        
     
        hold on;
        if vector_on_the_ring(2)*next_vector(1) - next_vector(2)*vector_on_the_ring(1) ~= 0 
            intersection2 = (next_vector(1)*value - vector_on_the_ring(1)*next_value)/(vector_on_the_ring(2)*next_vector(1) - next_vector(2)*vector_on_the_ring(1));
            if vector_on_the_ring(1) ~= 0
                intersection1 = (value - intersection2*vector_on_the_ring(2))/vector_on_the_ring(1);
                line([point(1), intersection1], [point(2), intersection2], 'Color', 'blue');
                line([next_point(1), intersection1], [next_point(2), intersection2], 'Color', 'blue');         
            elseif next_vector(1) ~= 0
                intersection1 = (next_value - intersection2*next_vector(2))/next_vector(1);
                line([point(1), intersection1], [point(2), intersection2], 'Color', 'blue');
                line([next_point(1), intersection1], [next_point(2), intersection2], 'Color', 'blue');
            else
                line([ point(1), next_point(1) ], [ point(2), next_point(2) ], 'Color', 'blue');
            end
        else
            line([ point(1), next_point(1) ], [ point(2), next_point(2) ], 'Color', 'blue');
        end
        %i
        %vector_on_the_ring
        %next_vector
        point
        next_point
        line([point(1), next_point(1)], [point(2), next_point(2)], 'Color', 'red');
        axis([-8 8 -8 8]);
    end
    hold off;
end




