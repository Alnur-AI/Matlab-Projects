% preparations
close all
clear
clc

% input data
A = [1 -10;
     -5 1];
B = [1 0;
     0 1];
f = [0;
     0];
p = [4;
     4];   
a = 1;
t_0 = 0;
x_0 = [-5;
       -5];
r_0 = 1;
coef_alpha = 2;
coef_s = 12;
time = 5;

num_of_iterations = 30;

%array for circle X_0
alpha_round = linspace(0, 2*pi, 500);

% plot X_0
figure(1);
plot(sin(alpha_round)*r_0+x_0(1), cos(alpha_round)*r_0+x_0(2), '.m');
hold on
axis([ -15, 15, -15, 15]);

% plot X_1
x1_lin = linspace(-coef_s/coef_alpha+1, coef_s/coef_alpha+1, 100);
x2_lin = linspace(1-sqrt(coef_s), 1+sqrt(coef_s), 100);
plot(x1_lin, ones(size(x1_lin))*(-sqrt(coef_s)+1), '.-g');
plot(x1_lin, ones(size(x1_lin))*(sqrt(coef_s)+1), '.-g');
plot(repmat([-coef_s/coef_alpha+1], size(x2_lin)), x2_lin, '.-g');
plot(ones(size(x2_lin))*(coef_s/coef_alpha+1), x2_lin, '.-g');

%plot P
figure(4)
hold on
x_lin = linspace(-1, 1, 500);
plot(cos(linspace(0, pi/2, 500))+a+p, sin(linspace(0, pi/2, 500))+a+p, '.g');
plot(cos(linspace(pi/2, pi, 500))-a+p, sin(linspace(pi/2, pi, 500))+a+p, '.g');
plot(cos(linspace(3*pi/2, 2*pi, 500))+a+p, sin(linspace(3*pi/2, 2*pi, 500))-a+p, '.g');
plot(cos(linspace(pi, 3*pi/2, 500))-a+p, sin(linspace(pi, 3*pi/2, 500))-a+p, '.g');
plot(x_lin*a+p, a+p+1, '.g')
plot(x_lin*a+p, -a+p-1, '.g')
plot(-a+p-1, x_lin*a+p, '.g')
plot(a+p+1, x_lin*a+p, '.g')

% starting enumeration of psi_0
alpha = linspace(0, 2*pi, num_of_iterations)';
psi_0 = [cos(alpha) sin(alpha)]';

% getting vector of couples [psi(t)_1, psi(t)_2, ...]
% [ x1 x2 ...
% [ y1 y2 ...
psi = @(t) expm(-A'.*t)*expm(A'.*t_0)*psi_0;

% finding start position
x_t_0 = x_0 + r_0*psi(t_0);

% finding optimal control for each psi(t)_i
l = @(t) B'*psi(t);
beta = @(t) acos( ([1 0]*l(t))./(([1 0]*l(t)).^2 + ([0 1]*l(t)).^2).^(1/2) );
%u = @(t) l(t) ./ ((([1 0]*l(t)).^2 + ([0 1]*l(t)).^2).^(1/2)) .* (a*(abs([1 0]*l(t)) + abs([0 1]*l(t))) ./ ((([1 0]*l(t)).^2 + ([0 1]*l(t)).^2).^(1/2)) + 1) + p;
num_points_u = 1;
if (rank(B) == 0)
    u = @(t) t.*zeros(size(psi_0)) + a + p + [0;1];
elseif (rank(B) == 2)
    u = @(t) a*sign([cos(beta(t)); sin(beta(t))]) + p + [cos(beta(t)); sin(beta(t))];
elseif (rank(B) == 1)
    u_2_matrix = a*sign([cos(beta(1)); sin(beta(1))]) + p + [cos(beta(1)); sin(beta(1))];
    u_2_matrix = [u_2_matrix(:,min(find(u_2_matrix(1,:) == min(u_2_matrix)))) u_2_matrix(:,min(find(u_2_matrix(1,:) == max(u_2_matrix))))];
    num_points_u = 2;
    %u_matrix = [linspace(u_2_matrix(1, 1), u_2_matrix(1, 2), num_points_u); linspace(u_2_matrix(2, 1), u_2_matrix(2, 2), num_points_u)]
    u_matrix = u_2_matrix;
    u = @(t) t.*zeros(size(u_matrix)) + u_matrix;
end
    
% finding the trajectory that hit the X_1
options = odeset('AbsTol', 1e-7, 'RelTol', 1e-5, 'Events',@(t,x)reached_the_X_1(t,x,coef_s, coef_alpha));%, 'MaxStep', 1e-6)
tspan = [t_0 t_0+time];
t_1_inf = t_0+time;
plot_x = [];

size_psi = 0;

% using ode45 for finding the trajectory that hit the X_1
for i = 1:size(psi_0, 2)
    for d = 1:num_points_u
        matrix = zeros(size(psi_0, 2), 1);
        matrix(i, 1) = 1;
        
        osobui_matrix_u = [];
        if (num_points_u > 1)
            osobui_matrix_u = zeros(size(u_matrix, 2), 1);
            osobui_matrix_u(d, 1) = 1;
        end

        %[t, x, te, xe, ie] = ode45( @(t, x) A*[x(1); x(2)] + B*u(t)*matrix + f, tspan, x_t_0*matrix, options);
        [t, x, te, xe, ie] = ode45( @(t,x) odefcn(t, x, A, B, f, @(t) u(t), matrix, osobui_matrix_u), tspan, x_t_0*matrix, options);
        
        size(x)
        
        %if (mod(i,2) < 1)

            if ~(max(( (coef_alpha*abs(x(1:end-1, 1)-2)<=coef_s) & ((x(1:end-1, 2)-2).^2<=coef_s) )))
                %plot each trajectory 
                figure(1);
                plot(x(:, 1), x(:, 2), '-c');
                hold on
            end

            figure(2);
            plot(t, x(:, 1), '-c');
            hold on

            figure(3);
            plot(t, x(:, 2), '-c');
            hold on
            
            if (num_points_u > 1)
                plot_y_u = [];
                for count = 1:length(t)
                    plot_y_u = [plot_y_u u(t(count))*osobui_matrix_u];
                end
            else
                plot_y_u = [];
                for count = 1:length(t)
                    plot_y_u = [plot_y_u u(t(count))*matrix];
                end
            end
            
            figure(4);
            plot_y_u = plot_y_u;
            plot(plot_y_u(1, :), plot_y_u(2, :), '*m');
            hold on

            figure(5);
            plot_y_u = plot_y_u;
            plot(t, plot_y_u(1, :), 'c');
            hold on

            figure(6);
            plot_y_u = plot_y_u;
            plot(t, plot_y_u(2, :), 'c');
            hold on

            plot_y_psi = [];
            for count = 1:length(t)
                plot_y_psi = [plot_y_psi psi(t(count))*matrix];
            end

            if ~(max(( (coef_alpha*abs(x(1:end-1, 1)-1)<=coef_s) & ((x(1:end-1, 2)-1).^2<=coef_s) )))
                figure(7);
                plot_y_psi = plot_y_psi;
                plot(plot_y_psi(1, :), plot_y_psi(2, :), 'c');
                hold on
                size_psi = size(plot_y_psi);
            else
                %figure(7);
                %plot_y_psi_1 = plot_y_psi(:,1:size_psi(2));
                %plot(plot_y_psi_1(1, :), plot_y_psi_1(2, :), 'c');
                %hold on
                figure(7);
                plot_y_psi = plot_y_psi;
                plot(plot_y_psi(1, :), plot_y_psi(2, :), 'c');
                hold on
                size_psi = size(plot_y_psi);
            end

            figure(8);
            plot_y_psi = plot_y_psi;
            plot(t, plot_y_psi(1, :), 'c');
            hold on

            figure(9);
            plot_y_psi = plot_y_psi;
            plot(t, plot_y_psi(2, :), 'c');
            hold on
        %end

        % store the optimal ones
        if((~isempty(te)) && (te < t_1_inf) && (~isnan(te)) && (~isempty(xe)) && (~isnan(xe(1))) && (~isnan(xe(2))))
            t_1_inf = te;
            x_res_inf = xe;
            plot_t = t;
            plot_x = x;
            x_00000 = x_t_0*matrix;
            
            l_for_T_r = psi(t_1_inf)*matrix/norm(psi(t_1_inf)*matrix)*(-1);
            norm_T_r = [0 0];
            if (abs(l_for_T_r(1)) > 0.3)
                norm_T_r(1) = sign(l_for_T_r(1));
            end
            if (abs(l_for_T_r(2)) > 0.3)
                norm_T_r(2) = sign(l_for_T_r(2));
            end
            T_r = cos(dot(l_for_T_r, norm_T_r));
            %T_r = l_for_T_r(1)*x_res_inf(1) + l_for_T_r(2)*x_res_inf(2);
            %T_r = abs(  T_r - ( l_for_T_r(1) + l_for_T_r(2) + abs(l_for_T_r(1))*coef_s/coef_alpha + abs(l_for_T_r(2))*sqrt(coef_s) )  );
        end
    end
end

figure(1);
xlabel('$x_1$','Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');

figure(2);
xlabel('$t$','Interpreter','latex');
ylabel('$x_1$','Interpreter','latex');

figure(3);
xlabel('$t$','Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');

figure(4);
xlabel('$u_1$','Interpreter','latex');
ylabel('$u_2$','Interpreter','latex');

figure(5);
xlabel('$t$','Interpreter','latex');
ylabel('$u_1$','Interpreter','latex');

figure(6);
xlabel('$t$','Interpreter','latex');
ylabel('$u_2$','Interpreter','latex');

figure(7);
xlabel('$\psi_1$','Interpreter','latex');
ylabel('$\psi_2$','Interpreter','latex');

figure(8);
xlabel('$t$','Interpreter','latex');
ylabel('$\psi_1$','Interpreter','latex');

figure(9);
xlabel('$t$','Interpreter','latex');
ylabel('$\psi_2$','Interpreter','latex');

if max(size(plot_x)) < 1
    plot_x = zeros(length(t), length(t));
    display('There is no solution.');
    return;
end

% plot optimal trajectory, u, psi
figure(1);
plot(plot_x(:, 1), plot_x(:, 2), '-k');
hold on

figure(2);
plot(plot_t, plot_x(:, 1), '-k');
hold on
axis([ min(plot_t) max(plot_t) min(plot_x(:, 1)) max(plot_x(:, 1))]);

figure(3);
plot(plot_t, plot_x(:, 2), '-k');
hold on
axis([ min(plot_t) max(plot_t) min(plot_x(:, 2)) max(plot_x(:, 2))]);

if (num_points_u > 1)
    plot_u = [];
    for count = 1:length(plot_t)
        plot_u = [plot_u u(plot_t(count))*osobui_matrix_u];
    end
else
    plot_u = [];
    for count = 1:length(plot_t)
        plot_u = [plot_u u(plot_t(count))*matrix];
    end
end


figure(4);
plot_u = plot_u;
plot(plot_u(1, :), plot_u(2, :), '*k');
hold on
axis equal

figure(5);
plot_u = plot_u;
plot(plot_t, plot_u(1, :), '-k');
hold on
xlim([0 1]);
xlim([0 max(plot_t)]);

figure(6);
plot_u = plot_u;
plot(plot_t, plot_u(2, :), '-k');
hold on
xlim([0 1]);
xlim([0 max(plot_t)]);

plot_psi = [];
for count = 1:length(plot_t)
    plot_psi = [plot_psi psi(plot_t(count))*matrix];
end

figure(7);
plot_psi = plot_psi;
plot(plot_psi(1, :), plot_psi(2, :), '-k');
hold on
xlim([-2 2]);
ylim([-2 2]);
%xlim([min(plot_psi(1,:)) max(plot_psi(1,:))])
%ylim([min(plot_psi(2,:)) max(plot_psi(2,:))])

figure(8);
plot_psi = plot_psi;
plot(plot_t, plot_psi(1, :), '-k');
hold on
xlim([0 1]);
xlim([0 max(plot_t)])

figure(9);
plot_psi = plot_psi;
plot(plot_t, plot_psi(2, :), '-k');
hold on
xlim([0 1]);
xlim([0 max(plot_t)])

figure(1);
plot(linspace(x_res_inf(1), x_res_inf(1) + l_for_T_r(1), 30), linspace(x_res_inf(2), x_res_inf(2) + l_for_T_r(2), 30), '.r')
hold on
plot(x_res_inf(1), x_res_inf(2), '.r');
plot(linspace(x_res_inf(1), x_res_inf(1) + norm_T_r(1), 30), linspace(x_res_inf(2), x_res_inf(2) + norm_T_r(2), 30), '.b')
plot(x_res_inf(1), x_res_inf(2), '.b');

display('I am done!');
fprintf('The minimum time found is equal to %s.\n', num2str(t_1_inf));
fprintf('The T_r is equal to %s.\n', num2str(T_r));

function [value, isterminal, direction] = reached_the_X_1(t, x, coef_s, coef_alpha)
    value = ~( (coef_alpha*abs(x(1)-1)<=coef_s) & ((x(2)-1).^2<=coef_s) ); %value of function
    isterminal = 1;   %end the function at the zero point
    direction = 0;  %all the zeros 
end

function dydt = odefcn(t, x, A, B, f, u, matrix, osobui_matrix_u)
dydt = zeros(2,1);
if (max(size(osobui_matrix_u)) == 0)
    dydt = A*[x(1); x(2)] + B*u(t)*matrix + f;
else
    dydt = A*[x(1); x(2)] + B*u(t)*osobui_matrix_u + f;
end
end