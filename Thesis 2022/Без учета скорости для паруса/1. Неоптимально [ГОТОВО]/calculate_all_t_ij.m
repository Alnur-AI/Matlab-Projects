%% clear everything
clc 
clear

n = 9;% количество уровней
m = 13;% количество точек на уровне

% Делим поверхность на сетку
Vx = linspace(-1000, 1000, m);
Vy = linspace(0, 2000, n);

% Параметры лодки и ветра
boat_k = 5;
boat_k1 = 10;
boat_vb = 5;
boat_m = 1000;
boat_S = 10;

% Начальное положение, скорость и начальное время
boat_x = 0;
boat_y = 0;
boat_pos = 0;
boat_speed = 0;
time = 0;


%% calculate V(i,j)
clc


%Заполним нашу функцию цены
V = zeros(n,m);
%V(n,:) = floor(rand(m,1)*10000);

for i = flip(1:(n-1))
    
    
    for j = 1:m
        
        minV = inf;
        
        fprintf('\n')
        fprintf('\n')
        fprintf('i j l time\n')
        fprintf('\n')
        for l = 1:m
            
            
            if( abs( Vx(j) - Vx(l) ) > 0.1 )
           
                [boat_x, boat_y, boat_pos , boat_speed, time] = transfer_between_dots(Vx(j),Vy(i),Vx(l),Vy(i+1),boat_k,boat_k1,boat_vb,boat_m,boat_S, boat_pos, boat_speed,time);

                if (minV > V(i+1,l) + time && l ~= j) 
                    minV = V(i+1,l) + time;
                end

                fprintf('%i %i %i %d\n', i,j,l,time)

                boat_pos = 0;
                boat_speed = 0;
                time = 0;
                
            else
                fprintf('%i %i %i %d\n', i,j,l,inf)
            end
        end
        V(i,j) = minV;
    end
    
end
%% calculate gamma(t)

gamma_time = zeros(2*m-1,1);
gamma_value = zeros(2*m-1,1);

boat_pos = 0;
    boat_speed = 0;

% left
for i = 1:m-1

    init_x = Vx(m);
    init_y = Vy(1);
    dest_x = Vx(i);
    dest_y = Vy(2);
    [boat_x, boat_y, boat_pos , boat_speed, gamma_time(i)] = transfer_between_dots(init_x,init_y, dest_x,dest_y, boat_k,boat_k1,boat_vb,boat_m,boat_S, boat_pos, boat_speed,time);
    gamma_value(i) = pi+ acot( (init_x-dest_x)/(init_y-dest_y) );

    boat_pos = 0;
    boat_speed = 0;
end

gamma_time(m) = inf;
gamma_value(m) = pi/2;

% right
for i = (1:m-1)

    init_x = Vx(1);
    init_y = Vy(1);
    dest_x = Vx(i+1);
    dest_y = Vy(2);
    [boat_x, boat_y, boat_pos , boat_speed, gamma_time(m+i)] = transfer_between_dots(init_x,init_y, dest_x,dest_y, boat_k,boat_k1,boat_vb,boat_m,boat_S, boat_pos, boat_speed,time);
    gamma_value(m+i)= acot( (init_x-dest_x)/(init_y-dest_y) );
    
    boat_pos = 0;
    boat_speed = 0;

end
gamma_value_angle = gamma_value*180/pi;
plot(gamma_value_angle, gamma_time, '*')

%% Find Optimal
%V = zeros(n,m);
j_opt = zeros(n,1);
pos_opt = zeros(n,1);
speed_opt = zeros(n,1);


j_opt(1) = floor(m/2);

for i = 1:n-1
    
    minV = inf;
    minV_ind = inf;
    
    fprintf('\n')
    fprintf('\n')
    fprintf('i j time V(i,j)\n')
    fprintf('\n')
    
    
    
    for j = 1:m
        
        if( abs( Vx(j) - Vx(j_opt(i)) ) < 0.1 )
            
            fprintf('%i %i  %d %d\n', i,j,inf, V(i+1,j))
            continue;
        end
        [boat_x, boat_y, boat_pos, boat_speed, time] = transfer_between_dots(Vx(j_opt(i)),Vy(i),Vx(j),Vy(i+1),boat_k,boat_k1,boat_vb,boat_m,boat_S, pos_opt(i), speed_opt(i),time);
        
        
        if( minV > time + V(i+1,j) && j ~= j_opt(i) )
            minV = V(i+1,j) + time;
            minV_ind = j;
            pos_opt(i+1) = boat_pos;
            speed_opt(i+1) = boat_speed;
        end
        
        
        fprintf('%i %i  %d %d\n', i,j,time, V(i+1,j))
        
        boat_pos = 0;
        boat_speed = 0;
        time = 0;
    end
    fprintf('Min time: %d    Min V(i,j): %d   Opt j: %i \n',minV - V(i+1,minV_ind) , minV, minV_ind)
    j_opt(i+1) = minV_ind;
end

%% plot 

plot( Vx(j_opt), Vy(1:n) )


