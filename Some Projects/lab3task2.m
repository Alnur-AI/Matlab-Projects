fun = @(x) cos(x) - x ./ pi;
grid = -5 : 0.1 : 5;
hold on;
plot(grid, cos(grid), 'b');
plot(grid, grid ./ pi, 'r');
plot(grid, grid ./ pi, 'r');
x0 = ginput(15);
for i = 1:size(x0, 1)
    plot(x0(i, 1), x0(i, 2), 'p', 'MarkerIndices', 1, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'red', 'MarkerSize', 8.5);
    approx_root = fzero(fun, x0(i, 1));
    discrepancy = fun(approx_root);
   
    
    fprintf('Если начальное приближение равно ')
    fprintf(num2str(x0(i, 1)));
    fprintf(', то приближенный корень уравнения равен ');
    fprintf(num2str(approx_root));
    fprintf('; невязка равна ');
    fprintf(num2str(discrepancy));
    fprintf('.\n');
    
    
    plot(approx_root, cos(approx_root), '-s', 'MarkerIndices', 1, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'green', 'MarkerSize', 8.5);
end
fprintf('\n');
hold off;
