%% Блок 1
f = @(x, y, c) 4 * sin(2 * c * pi * x) .* cos(c * 1.5 * pi * y) .* (1 - x .^ 2) .* y .* (1 - y) + c;

h = 0.05;
x = -1 : h : 1;
y = -1 : h : 1;
[X, Y] = meshgrid(x, y);

h_c = 0.5;
values_c = -1.5 : h_c : 1;

%% Блок 2
nFrames = size(values_c, 2);
mov(1:nFrames) = struct('cdata', [], 'colormap', []);
for i = 1:nFrames
    Z = f(X, Y, values_c(i));
    localmin_logic_matrix = findmax(-f(X, Y, values_c(i)));
    localmax_logic_matrix = findmax(f(X, Y, values_c(i)));
    
    surf(X, Y, f(X, Y, values_c(i)));
    axis([-1 1 0 1 -3 2]);
    hold on;
    Xmin = X .* localmin_logic_matrix;
    Ymin = Y .* localmin_logic_matrix;
    Zmin = Z .* localmin_logic_matrix;
    plot3(Xmin, Ymin, Zmin, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'red');
    Xmax = X .* localmax_logic_matrix;
    Ymax = Y .* localmax_logic_matrix;
    Zmax = Z .* localmax_logic_matrix;
    plot3(Xmax, Ymax, Zmax, 's', 'MarkerSize', 10, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'green');
    mov(i) = getframe();
    hold off;
    pause(0.5);
end

%% Блок 3

movie(mov);

%% Блок 4

c_fix = values_c(2);
level = 1.2;
Z = f(X, Y, c_fix);
contour(X, Y, Z, level);

%% Блок 5

v = VideoWriter('filevideo1.avi');
v.FrameRate = 0.5;
open(v);

for i = 1 : nFrames
    writeVideo(v, mov(i));
end
close(v);
save('filevideo1.mat', 'mov');
load('filevideo1.mat', '-mat');


function mx = findmax(A)
    n = size(A, 1);
    m = size(A, 2);
    
    Al = A;
    Al(:, 1) = [];
    Al(:, end + 1) = -Inf;
    
    Ar = A;
    Ar(:, end) = [];
    Ar = [ones(n, 1) * (-Inf), Ar];
    
    Au = A;
    Au(1, :) = [];
    Au(end + 1, :) = -Inf;
    
    Ad = A;
    Ad(end, :) = [];
    Ad = [ones(1, m) * (-Inf); Ad];
    
    Alu = Al;
    Alu(1, :) = [];
    Alu(end + 1, :) = -Inf;
    
    Ald = Al;
    Ald(end, :) = [];
    Ald = [ones(1, m) * (-Inf); Ald];
    
    Aru = Ar;
    Aru(1, :) = [];
    Aru(end + 1, :) = -Inf;
    
    Ard = Ar;
    Ard(end, :) = [];
    Ard = [ones(1, m) * (-Inf); Ard];
    
    mx = (A > Al) .* (A > Ar) .* (A > Au) .* (A > Ad);
    mx = mx .* (A > Alu) .* (A > Aru) .* (A > Ald) .* (A > Ard);
    
    if (A(1, 1) <= A(1, 2) || A(1, 1) <= A(2, 1) || A(1, 1) <= A(2,2))
        t = mx(:);
        t(t == 0) = inf;
        mx = reshape(t(:), n, m);
    end
end

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    





