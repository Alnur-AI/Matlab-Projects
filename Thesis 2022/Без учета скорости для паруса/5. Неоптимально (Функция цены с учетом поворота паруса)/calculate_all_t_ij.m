%% clear everything
clc
clear

n = 11;% количество уровней
m = 11;% количество точек на уровне
q = 10;% количество секторов (разделение строки j на q секторов)

% Делим поверхность на сетку
Vx = linspace(-1000, 1000, m);
Vy = linspace(0, 10000, n);

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

%% calculate gamma(t)

gamma_time = zeros(2*m-1,1);
gamma_value = zeros(2*m-1,1);

%boat_pos = 0;
%boat_speed = 0;

% left
tic
for i = 1:m-1

    init_x = Vx(m);
    init_y = Vy(1);
    dest_x = Vx(i);
    dest_y = Vy(2);
    [boat_x, boat_y, ~ , ~, gamma_time(i)] = transfer_between_dots(init_x,init_y, dest_x,dest_y, boat_k,boat_k1,boat_vb,boat_m,boat_S, boat_pos, boat_speed,time);
    gamma_value(i) = pi + acot( (init_x-dest_x)/(init_y-dest_y) );

    
    disp(i)
end

% middle
gamma_time(m) = inf;
gamma_value(m) = pi/2;

disp(m)

% right
for i = (1:m-1)

    init_x = Vx(1);
    init_y = Vy(1);
    dest_x = Vx(i+1);
    dest_y = Vy(2);
    [boat_x, boat_y, ~ , ~, gamma_time(m+i)] = transfer_between_dots(init_x,init_y, dest_x,dest_y, boat_k,boat_k1,boat_vb,boat_m,boat_S, boat_pos, boat_speed,time);
    gamma_value(m+i)= acot( (init_x-dest_x)/(init_y-dest_y) );
    
    
    disp(m+i)

end
toc
gamma_value_angle = gamma_value*180/pi;
%plot(gamma_value_angle, gamma_time, '*')

%% calculate V(i,j,p)
beta = 1;
V = zeros( n, m, q );
randvec = floor(rand(1,m)*100);
% phi(j) - случайные числа
%for p = 1:q
    %V(n,:,p) = randvec(:);
%end
% phi(j) - лучший вариант только середина
V(n,:,:) = 999;
V(n,floor(m/2),:) = 0;

minV_j = zeros(q,1);
Vijp_angle = zeros(q,1);

for i = flip(1:(n-1))
    for j = 1:m
        for p = 1:q
            
            minV_p = inf;
            minV_j(:) = inf;
            %Vijp_angle(:) = inf;
            
            
            for pf = 1:q
                
                minV_j(pf) = inf;
                
                for jf = 1:m
                    if( j ~= jf )
                        time = gamma_time( m-j+jf );
                        if (minV_j(pf) > V(i+1,jf,pf) + time) 
                            minV_j(pf) = V(i+1,jf,pf) + time;
                        end
                        %fprintf('%i %i %i %i %i %d %d %d\n', i,j,p,pf,jf,time,minV_j(pf),minV_p)
                    else
                        %fprintf('%i %i %i %i %i %d %d %d\n', i,j,p,pf,jf,inf,minV_j(pf),minV_p)
                    end
                end
                
                Vijp_angle(pf) = beta*abs(p-pf);
                
                %if (minV_p > beta*abs(p-pf) + minV_j(pf)) 
                %    minV_p = beta*abs(p-pf) + minV_j(pf);
                %end 
                
                
            end
            
            minV_p = min( Vijp_angle + minV_j );
            V(i,j,p) = minV_p;
            
            
            
        end
    end
end
   
%% Find Optimal V(i,j,p)

clc
j_opt = zeros(n,1);
p_opt = zeros(n,1);

j_opt(1) = floor(m/2);
p_opt(1) = floor( j_opt(1)*q/(m-1) );

minV_j = zeros(q,1);
minV_j_ind = zeros(q,1);

total_time = 0; % Общее время для плавания

for i = 1:n-1
    minV_j(:) = inf;
    minV_j_ind(:) = inf;
    minV_p = inf;
    minV_p_ind = inf;

    for pf = 1:q
        
        minV_j(pf) = inf;
        minV_j_ind(pf) = inf;
        
        for jf = 1:m
            
            time = gamma_time(m-j_opt(i)+jf);
            if(minV_j(pf) > time + V(i+1,jf,pf) && jf ~= j_opt(i))
                minV_j(pf) = time + V(i+1,jf,pf);
                minV_j_ind(pf) = jf;
            end
        end
        
        
        if(minV_p > beta*abs(floor( j_opt(i)*q/(m-1) ) - pf) + minV_j(pf))
            minV_p = beta*abs(floor( j_opt(i)*q/(m-1) ) - pf) + minV_j(pf);
            minV_p_ind = pf;
        end
        %{
        if(minV_p > beta*abs(p_opt(i) - pf) + minV_j(pf))
            minV_p = beta*abs(p_opt(i) - pf) + minV_j(pf);
            minV_p_ind = pf;
        end
        %}
        %disp(minV_j_ind)
    end
    
    total_time = total_time + gamma_time(m-j_opt(i)+minV_j_ind(minV_p_ind));
    
    j_opt(i+1) = minV_j_ind((minV_p_ind));
    p_opt(i) = (minV_p_ind);
        
end
%% Find Optimal V(i,j,p) without stoping ship

clc
j_opt = zeros(n,1);
p_opt = zeros(n,1);

j_opt(1) = floor(m/2);
p_opt(1) = floor( j_opt(1)*q/(m-1) );

minV_j = zeros(q,1);
minV_j_ind = zeros(q,1);

pos_opt = zeros(q,n);
speed_opt = zeros(q,n);

total_time = 0;
min_time = zeros(q,1);

for i = 1:n-1
    minV_j(:) = inf;
    minV_j_ind(:) = inf;
    minV_p = inf;
    minV_p_ind = inf;

    for pf = 1:q
        
        minV_j(pf) = inf;
        minV_j_ind(pf) = inf;
        
        for jf = 1:m
            
            %time = gamma_time(m-j_opt(i)+jf);
            [boat_x, boat_y, boat_pos, boat_speed, time] = transfer_between_dots(Vx(j_opt(i)),Vy(i),Vx(j),Vy(i+1),boat_k,boat_k1,boat_vb,boat_m,boat_S,  pos_opt(pf,i), speed_opt(pf,i),time);
            
            if(minV_j(pf) > time + V(i+1,jf,pf) && jf ~= j_opt(i))
                minV_j(pf) = time + V(i+1,jf,pf);
                minV_j_ind(pf) = jf;
                
                pos_opt(pf,i+1) = boat_pos;
                speed_opt(pf,i+1) = boat_speed;
                min_time(pf) = time;
            end
            %boat_pos = 0;
            %boat_speed = 0;
            time = 0;
        end
        
        
        if(minV_p > beta*abs(floor( j_opt(i)*q/(m-1) ) - pf) + minV_j(pf))
            minV_p = beta*abs(floor( j_opt(i)*q/(m-1) ) - pf) + minV_j(pf);
            minV_p_ind = pf;
        end
        
    end
    
    %total_time = total_time + gamma_time(m-j_opt(i)+minV_j_ind(minV_p_ind));
    total_time = total_time + min_time(minV_p_ind);
    
    j_opt(i+1) = minV_j_ind((minV_p_ind));
    p_opt(i) = (minV_p_ind);
        
end
%% plot V(i,j,p)

figure 
plot(gamma_value_angle, gamma_time, '*')
figure 
plot( Vx(j_opt), Vy(1:n) )
%% calculate V(i,j)
clc


%Заполним нашу функцию цены
V = zeros(n,m);
V(n,:) = floor(rand(m,1)*100);
V(n,:) = 999;
V(n,floor(m/2)+1) = 0;

for i = flip(1:(n-1))
    
    
    for j = 1:m
        
        minV = inf;
        
        fprintf('\n')
        fprintf('\n')
        fprintf('i j l time\n')
        fprintf('\n')
        for l = 1:m
            
            
            if( abs( Vx(j) - Vx(l) ) > 0.1 )
           
                %[boat_x, boat_y, boat_pos , boat_speed, time] = transfer_between_dots(Vx(j),Vy(i),Vx(l),Vy(i+1),boat_k,boat_k1,boat_vb,boat_m,boat_S, boat_pos, boat_speed,time);
                time = gamma_time( m-j+l );
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

%% Find Optimal 
%V = zeros(n,m);
clc
j_opt = zeros(n,1);
pos_opt = zeros(n,1);
speed_opt = zeros(n,1);
total_time = 0; % Общее время для плавания

j_opt(1) = floor(m/2)+1;

for i = 1:n-1
    
    minV = inf;
    minV_ind = inf;
    
    fprintf('\n')
    fprintf('\n')
    fprintf('i j time V(i,j) overall\n')
    fprintf('\n')
    if (mod(i,2)==1)
        jloop = 1:m;
    end
    if (mod(i,2)==0)
        jloop = flip(1:m);
    end
    for j = jloop
        
        % Расчет времени при моментальном торможении
        %time = gamma_time(m - j_opt(i) + j);
        
        if( abs( Vx(j) - Vx(j_opt(i)) ) < 0.1 )
            
            fprintf('%i %i  %d %d\n', i,j,time, V(i+1,j))
            continue;
        end
        % Расчет при несбрасывании скорости
        [boat_x, boat_y, boat_pos, boat_speed, time] = transfer_between_dots(Vx(j_opt(i)),Vy(i),Vx(j),Vy(i+1),boat_k,boat_k1,boat_vb,boat_m,boat_S,  pos_opt(i), speed_opt(i),time);
        
        
        if( minV > time + V(i+1,j) && j ~= j_opt(i) )
            minV = V(i+1,j) + time;
            minV_ind = j;
            pos_opt(i+1) = boat_pos;
            speed_opt(i+1) = boat_speed;
        end
        
        
        %fprintf('%i %i  %d %d %d\n', i,j,time, V(i+1,j), V(i+1,j) + time)
        fprintf('%i %i  %d\n', i,j, V(i+1,j) + time)
        
        %boat_pos = 0;
        %boat_speed = 0;
        time = 0;
    end
    %fprintf('Min time: %d  Min V(i,j): %d\nOpt j: %i  Overall: %d\n',minV - V(i+1,minV_ind) , V(i+1,minV_ind), minV_ind, minV)
    fprintf('Min time: %d  Min V(i,j): %d\nOpt j: %i  Overall: %d\n',minV - V(i+1,minV_ind) , V(i+1,minV_ind), minV_ind, minV)
    
    total_time = total_time + minV - V(i+1,minV_ind);
    j_opt(i+1) = minV_ind;
end

%% plot 


figure 
plot(gamma_value_angle(22:41), gamma_time(22:41), '*')
figure 
plot( Vx(j_opt), Vy(1:n) )
%plot( Vx(j_opt(1:floor(n/2))), Vy(1:floor(n/2)) )

