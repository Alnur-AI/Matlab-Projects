A = 3;
a = 1;
delta = pi/sqrt(2);
B = 2;
b = 0.57;

f = @(t) sin(t .* a + delta) .* A;
g = @(t) sin(t .* b) .* B;

%f = @(t) sin(t);
%g = @(t) cos(t);

t0 = 0;
t1 = 14;

N = 15;



%_________________________________________________



dots_on_the_line = getEqual(f, g, t0, t1, N)



%_________________________________________________
hold on;
h_uniform = (t1 - t0) / (N - 1);
uniform_setka = t0 : h_uniform : t1;
x = f(uniform_setka');
y = g(uniform_setka');

uniform_matrix(1 : N, 1) = x;
uniform_matrix(1 : N, 2) = y;

uniform_distance_matrix = pdist2(uniform_matrix, uniform_matrix);

sum = 0;
for j = 1 : N - 1
    sum = sum + uniform_distance_matrix(j, j + 1);
end

Average_distance_with_uniform_mesh = sum / (N - 1)

plot(x', y', 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'green', 'MarkerSize', 4);
axis equal
hold off;



%_________________________________________________



function [ array_of_dots ] = getEqual(f, g, t0, t1, N)
    little_segments_num = 5000; 
    h = (t1 - t0) / little_segments_num;
    setka = t0 : h : t1;
    
    matrix(1 : little_segments_num + 1, 1) = f(setka');
    matrix(1 : little_segments_num + 1, 2) = g(setka');
    distance_matrix = pdist2(matrix, matrix);

    dots_on_t0t1(1 : N - 1) = setka(1 : N - 1);
    dots_on_t0t1(N) = setka(end); 
    dots_on_t0t1_idxs(1 : N - 1) = 1 : N - 1;
    dots_on_t0t1_idxs(N) = size(setka, 2);
    
    for i = 1 : N -1
        array_of_distances(i) = distance_matrix(dots_on_t0t1_idxs(i), dots_on_t0t1_idxs(i + 1));
    end
    
    while max(array_of_distances) - min(array_of_distances) > 0.01
        idxs = find(array_of_distances == max(array_of_distances));
        idx = idxs(end);
        if idx ~= 1 
            dots_on_t0t1(idx) = dots_on_t0t1(idx) + h;
            dots_on_t0t1_idxs(idx) = dots_on_t0t1_idxs(idx) + 1;
            array_of_distances(idx) = distance_matrix(dots_on_t0t1_idxs(idx), dots_on_t0t1_idxs(idx + 1));
            array_of_distances(idx - 1) = distance_matrix(dots_on_t0t1_idxs(idx - 1), dots_on_t0t1_idxs(idx));
        else
            break;
        end
    end
    
    array_of_dots(1:N, 1) = f(dots_on_t0t1');
    array_of_dots(1:N, 2) = g(dots_on_t0t1');
    
    plot(f(setka), g(setka), '-o', 'MarkerIndices', dots_on_t0t1_idxs, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'red', 'MarkerSize', 8);
    hold on;
    h_uniform = (t1 - t0) / (N - 1);
    uniform_setka = t0 : h_uniform : t1;
    plot(f(uniform_setka), g(uniform_setka), 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'green', 'MarkerSize', 4);
    axis equal
    hold off;
    
    Required_distance = array_of_distances(1)
end








