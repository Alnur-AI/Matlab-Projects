points = [-1, -1; 5, 8; 13, 8; -4, 9; -2.5, 5; 0, -2.5];
V = 10;
L = 8;

viewPossible(points, V, L);

function [] = viewPossible(points, V, L)
    a = -5;
    b = 15;
    h = 0.25;
    [X, Y] = meshgrid(a : h : b, a : h : b);
    N = size(points, 1);
    Z = zeros(size(Y, 1), size(Y, 2));
    for i = 1:N
        Zi = ones(size(Y, 1), size(Y, 2)) .* V;
        station_coordinates = points(i, :);
        Zi = Zi ./ (1 + sqrt((Y - station_coordinates(2)).^2 + (X - station_coordinates(1)).^2));
        Z = Z + Zi;
    end
    M = contourf(X, Y, Z, [L, L])
    hold on;
    plot(points(:, 1), points(:, 2), 'p', 'MarkerIndices', 1:N, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'red', 'MarkerSize', 10);
    hold off;
    if size(M, 2) == M(2, 1) + 1
        title('односвязная');
    else 
        title('не односвязная');
    end
end









