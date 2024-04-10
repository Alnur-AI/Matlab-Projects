clc
clear

% Данная программа вычисляет универсально оптимальную траекторию
% Для уравнения x'= Ax + u, x - вектор второго размера
% На отрезке времени [t0, T]
% M_0 - Начальное множество
% U - Область управления
% Функционал: J = min (t - t0)

%% Часть 1: Задание оптимальных множеств, матрицы A и отрезка времени
% Легче всего задать замкнутое множество F через полярные координаты.
% Причём центр масс множества F лежит в (0,0)
% Если же это не так, центр масс добавляется в конце к каждой координате.
% Так же очевидно что 
% theta = linspace(0,2pi) - угол от 0 до 2pi
% r = r(theta) - радиус
% f0 = (x0,y0) - координаты центра масс 
% Получатся, что точка (x,y), лежащая на границе F можно вычислить
% следующим образом: 
% x = r*cos(theta) + x0;
% y = r*sin(theta) + y0;
% ПРИМЕЧАНИЕ: функции r(theta) заданы в отдельных функциях r_F,  

clc
clear

% Матрица
A = [1,1;-1,0];

% Временной промежуток
N_t = 70;%Разбиение на интервала t
t0 = 0;
T = 1;
t = linspace(t0,T,N_t);

% Вектор psi
N_psi = 70;% Разбиение на интервалы psi, M_0, U, M_1
theta = linspace(0,2*pi,N_psi);

% Начальное Множество
cM_0 = [1;0];% Координаты центра масс множества M_0

i = 1:N_psi;
r = zeros(N_psi,1);
r(i) = r_M_0(theta(i));

M_0 = zeros(N_psi,2);% Набор из точек что лежат на границе M_0
M_0(:,1) = r(:).*cos(theta(:)) + cM_0(1);
M_0(:,2) = r(:).*sin(theta(:)) + cM_0(2);

% Множество Управления
cU = [0;0];% Координаты центра масс множества U

i = 1:N_psi;
r = zeros(N_psi,1);
r(i) = r_U(theta(i));

U = zeros(N_psi,2);% Набор из точек что лежат на границе U
U(:,1) = r(:).*cos(theta(:)) + cU(1);
U(:,2) = r(:).*sin(theta(:)) + cU(2);

% Строим множества M_0 и U
figure (1)
hold on
plot(M_0(:,1),M_0(:,2),'-*');
plot(U(:,1),U(:,2),'-*');
legend('M_0','U')
hold off

clear i r

%% Часть 2: Вычисление опорной функции C(X(T),psi)
% Отсюда можно найти множество достижимости из следующих формул:
% CX(psi,t) = C(x(t),psi) = C(M_0, expm((t-t0)*A')*psi) + integral(   C(U,expm((t-s)*A')*psi )    ,t0,t,ds)

CX = zeros(N_psi,N_t);

% Для каждого отдельного вектора psi(theta)
for j = 1:N_psi
    psi_j = [cos(theta(j)); sin(theta(j))];
    
    % Для каждой переменной t(i)
    for i = 1:N_t
        
        % Вычисляем интеграл от t0 до t
        % integral(t0,t, C(U, expm((t-s)*A')*psi ),  ds );
        %%{
        CU = zeros(i,1);
        for k = 1:i
            CU(k) = C(U, expm((t(i)-t(k)).*transpose(A))*psi_j);
        end
        integral_CX = trapz(t(k),CU(k),2);
        %%}
        
        
        % Вычисляем интеграл от t0 до T
        % integral(t0,T, C(U, expm((t-s)*A')*psi ),  ds );
        %{
        CU = zeros(N_t,1);
        for k = 1:N_t
            CU(k) = C(U, expm((t(i)-t(k)).*transpose(A))*psi_j);
        end
        integral_CX = trapz(t,CU);
        %}
        
        
        CX(j,i) = C(M_0, expm((t(i)-t0).*transpose(A))*psi_j) + integral_CX;
        fprintf('%d %d\n', j , i);
    end

end

%% Часть 3.1: Строим график CX в декартовых координатах CX(theta,t)
% Так как мы нашли чему равна опорная функция C(x(t),psi)
% в предыдущей части, то можно построить трехмерный график
% по оси Ox - время, Oy - угол, Oz - значение опорной функции в этой точке
[XX,YY] = meshgrid(t,theta);

figure (1);
surf(XX,YY,CX)
xlabel('Time')
ylabel('psi(Angle)')
zlabel('CX')

clear XX YY

%% Часть 3.2: Строим график CX в цилиндрических координатах CX(theta,t)
% Не трудно заметить что функция CX(psi,t) имеет период 2pi по оси Ox
% Следовательно, представить эту функцию в 
% цилиндрических координатах может оказаться
% очень удобным
% В нашем случае:
% r = t
% angle = theta;
[XX,YY] = meshgrid(theta,t);
[xx,yy,zz] = pol2cart(XX,YY,CX);

figure (2)
surf(xx,yy,zz');
ylabel('x')
xlabel('y')
zlabel('CX(t*cos(l),t*sin(l))')

clear XX YY xx yy zz

%% Часть 4.1: Найдем x1(t) x2(t) по формуле:
%c = CX;% Тогда c_i = CX(psi_i,x);
% x1(i)*cos(theta(j))   + x2(i)*sin(theta(j)) = c(j,i)
% x1(i)*cos(theta(j+1)) + x2(i)*sin(theta(j+1)) = c(j+1,i)
% Таким образом из всевозможных траекторий (x1(t),x2(t))
% и всевозможных множеств достижимости от t0 до T
% находим следующую переменную:
% xfull(psi,t) в нашем случае она представляет собой 
% Трубу достижимости от t0 до T


xfull = zeros(N_psi,N_t,2);

for j = 1:N_psi-1
   for i = 1:N_t
       xfull(j,i,:) = [cos(theta(j)) sin(theta(j)) ; cos(theta(j+1)) sin(theta(j+1))  ]\[CX(j,i);CX(j+1,i)];
   end
end
i = 1:N_t;
xfull(N_psi,i,:) = xfull(1,i,:);

%% Часть 4.2: Вычислим M_1 - множество достижимости
M_1 = zeros(N_psi,2);

M_1(:,1) = xfull(:,end,1);
M_1(:,2) = xfull(:,end,2);

%% Часть 5.1: Строим в 2D всевозможные траектории (x1(t),x2(t)) от M_0 до M_1   
figure (3)
hold on
for i = 1:N_t
    plot( xfull(:,i,1), xfull(:,i,2) );
    pause(0.01)
end
hold off

figure (4)
hold on
    plot(xfull(:,1,1),xfull(:,1,2))
    pause(0.1)
    for j = 1:N_psi
        plot( xfull(j,:,1), xfull(j,:,2) );
        pause(0.01)
    end
    plot(xfull(:,end,1),xfull(:,end,2))
hold off

%% Часть 5.2: Строим в 3D всевозможные траектории (x1(t),x2(t),t) от M_0 до M_1  

figure (5)
view(3)

hold on
    plot3(M_0(:,1),M_0(:,2),ones(N_psi,1).*t(1));
    i = 1:N_t;
    plot3( xfull(:,i,1), xfull(:,i,2), ones(N_t,1).*t(i) );
hold off

figure (6)
view(3)
hold on

    plot3(M_0(:,1),M_0(:,2),ones(N_psi,1).*t(1));
    plot3( xfull(:,1,1), xfull(:,1,2), ones(N_t,1).*t(1) );
    
    for j = 1:N_psi
        plot3( xfull(j,:,1), xfull(j,:,2), ones(N_t,1).*t(:) );
    end
    plot3( xfull(:,N_t,1), xfull(:,N_t,2), ones(N_t,1).*t(N_t) );
hold off

%% Часть 6.1: Строим в 2D M_0 и M_1
figure (7)
hold on
    plot( M_0(:,1), M_0(:,2), '-*', mean(M_0(:,1)), mean(M_0(:,2)), '*' );
    plot( M_1(:,1), M_1(:,2), '-*', mean(M_1(:,1)), mean(M_1(:,2)), '*' );
    legend('M_0','Центр масс M_0','M_1','Центр масс M_1');
hold off





