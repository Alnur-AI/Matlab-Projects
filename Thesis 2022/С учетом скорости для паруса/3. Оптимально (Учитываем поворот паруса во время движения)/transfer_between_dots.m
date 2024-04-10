function [curr_x, curr_y, curr_pos, curr_speed, curr_time] = transfer_between_dots(init_x,init_y,dest_x,dest_y,boat_k,boat_k1,boat_vb,boat_m,boat_S, boat_pos,  boat_speed,time)

% Настройки расчетов
Tdelay = 1;%с (Временной промежуток для проверки паруса и прибытие до места назначения)

% INPUT DATA (А так же нужно отправить параметры корабля)
init_speed = boat_speed;%м/с (Начальная скорость  на прямой)

% Свойства Моря и ветра
k = boat_k;%        (Коэффициент полезного действия силы ветра)
k1 = boat_k1;%       (Коэффициент сопротивления воды)
vb = boat_vb;%м/с     (Скорость ветра)
av = pi/2;

% Судно
m = boat_m;%кг    (Масса судна)
S = boat_S;%м^2     (Площадь паруса)

% Временные переменные для расчета пути
curr_time = time;%                         Сколько времцени прошло с начала движения
curr_gamma = acot( (init_x-dest_x)/(init_y-dest_y) );% Текущий угол поворота судна. 
curr_pos = boat_pos;
curr_speed = init_speed;%                Текущая скорость, меняется каждую итерацию

curr_x = init_x;% Где корабль находился по х
curr_y = init_y;% Где корабль находился по у

% Пока не дойдем до точки достижения - будем считать
while( dest_y - curr_y > 0  )

    
    % Углы (С учетом скорости)
    angle_CBA = pi - (av - curr_gamma);
    
    % С корректировкой
    vb = sqrt(boat_vb^2 + curr_speed^2 -2*boat_vb*curr_speed*cos(angle_CBA) );
    delta = acos( (boat_vb^2 + vb^2-curr_speed^2)/(2*boat_vb*vb) );
    
    % Без корретировки
    %vb = boat_vb;
    %delta = 0;
    
    av =     90*(pi/180) + delta ;     % (Угол ветра относительно ОНСК)
    gamma = curr_gamma;    % (Угол направления движения судна относительно ОХ в ОНСК)
    theta = (av-gamma)/2 -delta/2;         % (Угол между парусом и носом судна)
    alpha =  theta;
    fprintf('vb = %f, av = %f, theta = %f\n',vb,av*180/pi,theta*180/pi);
    %disp((pi/2 - 2*theta - gamma+delta)*180/pi);
    
    % если curr_alpha > 90
    if (gamma < 0)
        av =     90*(pi/180) - delta;
        gamma = pi + gamma;
        theta = abs(gamma - av)/2;
        alpha = theta;
        disp((gamma)*180/pi);
    end
    
    
    
    %{
    % Углы (Параметры, без учета скорости)
    av = 90*(pi/180);
    gamma = curr_gamma;    % (Угол направления движения судна относительно ОХ в ОНСК)
    theta = abs(av-gamma)/2;
    alpha = abs( abs(av-gamma) - theta ); % (Угол между парусом и ветром)
    
    
    % если curr_alpha > 90
    if (gamma < 0)
        gamma = pi + gamma;
        theta = abs(gamma - av)/2;
        alpha = theta;
    end
    %}
    
    
    
    % Начальные условия 
    x10 = curr_pos;%м       (Начальное положение на прямой)
    x20 = curr_speed;%м/с     (Начальная скорость  на прямой)
    tspan  = [curr_time curr_time + Tdelay];%с (Промежуток времени)
    x0 = [x10 x20];

    % УРАВНЕНИЕ
    opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
    [t,z] = ode45(@(t,z) odefcn(t,z,k1,vb,alpha,theta,k,S,m), tspan, x0,opts);

    % Переход в общую систему координат (ОНСК)
    x0 = curr_x;
    y0 = curr_y;
    x = z(:,1)*cos(gamma)+x0;
    y = z(:,1)*sin(gamma)+y0;  

    % Фиксируем новые значения
    curr_x = x(end);
    curr_y = y(end);
    curr_pos = z(end,1);
    curr_speed = z(end,2);
    curr_time = curr_time + Tdelay;
    
    %if(gamma > pi/2)
        %fprintf ('distance: %d\n', dest_y - curr_y)
    %end
    %plot(z(:,1),z(:,2));

end

% OUTPUT DATA
alpha*180/pi;
theta*180/pi;
gamma*180/pi;
av*180/pi;

curr_x = x(end);
curr_y = y(end);
curr_pos = z(end,1);
curr_speed = z(end,2);
%curr_time;

end

