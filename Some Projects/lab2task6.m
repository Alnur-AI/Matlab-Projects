%% Блок 1
f = @(x) sin(x); % переменная-функция
h = 0.01; % шаг сетки
x = 0 : h : 27; % сетка

%% Блок 2
y = f(x);
idxmax_array = find(y == max(y)); % индекс элемента сетки, в котором достигается глобальный максимум функции
idxmax = idxmax_array(1);
localmin_logic_arr = islocalmin(y);

if (y(1) < y(2))
    localmin_logic_arr(1) = 1;
end

if (y(end) < y(end - 1))
    localmin_logic_arr(end) = 1;
end

localmin_indices = find(localmin_logic_arr ~= 0); % индексы элементов сетки, в которых достигаются локальные минимумы функции

plot(x, f(x), '-ob', 'MarkerIndices', localmin_indices, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'green', 'MarkerSize', 10);
hold on;
plot(x(idxmax), y(idxmax), '-p', 'MarkerIndices', 1, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'red', 'MarkerSize', 10);

minabs = min(abs(localmin_indices - idxmax(1)));
idx = find(abs(localmin_indices - idxmax(1)) == minabs);
min_nearest_to_max_idx = localmin_indices(idx(1));

if idxmax > min_nearest_to_max_idx & h > 0
    h = -h;
end
if idxmax < min_nearest_to_max_idx & h < 0
    h = -h;
end

setka = x(idxmax) : h : x(min_nearest_to_max_idx);
comet(setka, f(setka));
hold off;


