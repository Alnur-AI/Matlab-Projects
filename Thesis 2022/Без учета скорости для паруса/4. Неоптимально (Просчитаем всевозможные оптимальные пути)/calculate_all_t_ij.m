%% clear everything
clc
clear

n = 10;% ���������� �������
m = 100;% ���������� ����� �� ������

% ����� ����������� �� �����
Vx = linspace(-10000, 10000, m);
Vy = linspace(0, 5000, n);

% ��������� ����� � �����
boat_k = 10;
boat_k1 = 10;
boat_vb = 5;
boat_m = 1000;
boat_S = 10;

% ��������� ���������, �������� � ��������� �����
boat_x = 0;
boat_y = 0;
boat_pos = 0;
boat_speed = 0;
time = 0;

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
    gamma_value(i) = pi + acot( (init_x-dest_x)/(init_y-dest_y) );

    boat_pos = 0;
    boat_speed = 0;
    
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
    [boat_x, boat_y, boat_pos , boat_speed, gamma_time(m+i)] = transfer_between_dots(init_x,init_y, dest_x,dest_y, boat_k,boat_k1,boat_vb,boat_m,boat_S, boat_pos, boat_speed,time);
    gamma_value(m+i)= acot( (init_x-dest_x)/(init_y-dest_y) );
    
    boat_pos = 0;
    boat_speed = 0;
    
    disp(m+i)

end
gamma_value_angle = gamma_value*180/pi;
plot(gamma_value_angle, gamma_time, '*')

%% calculate V(i,j)
clc


%�������� ���� ������� ����
V = zeros(n,m);
V(n,:) = floor(rand(m,1)*100);
V(n,:) = 999;
V(n,floor(m/2)) = 0;

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

%% Find Optimal (by sorting)
sortV
j_opt = zeros(n,1);
total_time = 0; % ����� ����� ��� ��������

% ��������� �����
j_opt(1) = floor(m/2);

for i = 1:n-1
    
    minV = inf;
    minV_ind = inf;
    
    fprintf('\n')
    fprintf('\n')
    fprintf('i j time V(i,j)\n')
    fprintf('\n')
    for j = 1:m
        
        time = gamma_time(m - j_opt(i) + j);
        
        % ���� ����� �� ����� 90 ��������
        if( abs( Vx(j) - Vx(j_opt(i)) ) < 0.1 )
            fprintf('%i %i  %d %d\n', i,j,time, V(i+1,j))
            continue;
        end
   
        
        if( minV > time + V(i+1,j) && j ~= j_opt(i) )
            minV = V(i+1,j) + time;
            minV_ind = j;
        end
        
        
        fprintf('%i %i  %d %d\n', i,j,time, V(i+1,j))
        
        boat_pos = 0;
        boat_speed = 0;
        time = 0;
    end
    fprintf('Min time: %d    Min V(i,j): %d   Opt j: %i \n',minV - V(i+1,minV_ind) , V(i+1,minV_ind), minV_ind)
    total_time = total_time + minV - V(i+1,minV_ind);
    j_opt(i+1) = minV_ind;
end

%% Find Optimal 
%V = zeros(n,m);
clc
j_opt = zeros(n,1);
total_time = 0; % ����� ����� ��� ��������

j_opt(1) = floor(m/2);

for i = 1:n-1
    
    minV = inf;
    minV_ind = inf;
    
    fprintf('\n')
    fprintf('\n')
    fprintf('i j time V(i,j) overall\n')
    fprintf('\n')
    for j = 1:m
        
        time = gamma_time(m - j_opt(i) + j);
        
        if( abs( Vx(j) - Vx(j_opt(i)) ) < 0.1 )
            
            fprintf('%i %i  %d %d\n', i,j,time, V(i+1,j))
            continue;
        end
        %[boat_x, boat_y, boat_pos, boat_speed, time] = transfer_between_dots(Vx(j_opt(i)),Vy(i),Vx(j),Vy(i+1),boat_k,boat_k1,boat_vb,boat_m,boat_S, boat_pos, boat_speed,time);
        
        
        if( minV > time + V(i+1,j) && j ~= j_opt(i) )
            minV = V(i+1,j) + time;
            minV_ind = j;
        end
        
        
        %fprintf('%i %i  %d %d %d\n', i,j,time, V(i+1,j), V(i+1,j) + time)
        fprintf('%i %i  %d\n', i,j, V(i+1,j) + time)
        
        boat_pos = 0;
        boat_speed = 0;
        time = 0;
    end
    %fprintf('Min time: %d  Min V(i,j): %d\nOpt j: %i  Overall: %d\n',minV - V(i+1,minV_ind) , V(i+1,minV_ind), minV_ind, minV)
    fprintf('Min time: %d  Min V(i,j): %d\nOpt j: %i  Overall: %d\n',minV - V(i+1,minV_ind) , V(i+1,minV_ind), minV_ind, minV)
    
    total_time = total_time + minV - V(i+1,minV_ind);
    j_opt(i+1) = minV_ind;
end

%% plot 

figure 
surf(V)

figure 
plot(gamma_value_angle, gamma_time, '*')
figure 
plot( Vx(j_opt), Vy(1:n) )


