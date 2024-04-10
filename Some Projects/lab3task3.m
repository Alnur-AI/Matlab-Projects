f = @(x) fun(x);

a = 0.01;
n = 9001;
h = (2 * a) / n;

abscissa_axis = -a : h : a;
ordinate_axis = zeros(1, length(abscissa_axis));

%options = optimset('FunValCheck','on');

for i = 1:length(abscissa_axis)
    if abscissa_axis(i) == 0
        ordinate_axis(i) = 0;
    else
        ordinate_axis(i) = fzero(f, abscissa_axis(i));
        if isnan(ordinate_axis(i))
            if abscissa_axis(i) > 0
                ordinate_axis(i) = 1 / pi;
            else
                ordinate_axis(i) = -(1 / pi);
            end
        end
        if sign(ordinate_axis(i)) ~= sign(abscissa_axis(i))
            ordinate_axis(i) = ordinate_axis(i) * (-1);
        end
        num = round(1 / (pi * ordinate_axis(i)));
        
        if num + 1 ~= 0
            current_difference = abs(ordinate_axis(i) - abscissa_axis(i));
            while abs(   1 / ((num + 1) * pi) - abscissa_axis(i)   ) < current_difference 
                current_difference = abs(   1 / ((num + 1) * pi) - abscissa_axis(i)   );
                num = num + 1;
                if num == 0
                    break;
                end
            end
            ordinate_axis(i) = 1 / (num * pi);
        end
        if num - 1 ~= 0
            current_difference = abs(ordinate_axis(i) - abscissa_axis(i));
            while abs(   1 / ((num - 1) * pi) - abscissa_axis(i)   ) < current_difference 
                current_difference = abs(   1 / ((num - 1) * pi) - abscissa_axis(i)   );
                num = num - 1;
                if num == 0
                    break;
                end
            end
            ordinate_axis(i) = 1 / (num * pi);
        end
    end
end    

fzero(f, 1.124)

plot(abscissa_axis, ordinate_axis);

function [ function_value ] = fun1(x, dot_near_zero)
    if abs(x) <= abs(dot_near_zero) 
        function_value = dot_near_zero * sin(1 / dot_near_zero);
    else
        function_value = x .* sin(1 ./ x);
    end
end

function [ function_value ] = fun(x)
    if x == 0 
        function_value = 0;
    else
        function_value = x .* sin(1 ./ x);
    end
end


