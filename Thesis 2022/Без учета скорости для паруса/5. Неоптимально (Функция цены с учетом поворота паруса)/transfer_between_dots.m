function [curr_x, curr_y, curr_pos, curr_speed, curr_time] = transfer_between_dots(init_x,init_y,dest_x,dest_y,boat_k,boat_k1,boat_vb,boat_m,boat_S, boat_pos,  boat_speed,time)

% ��������� ��������
Tdelay = 1;%� (��������� ���������� ��� �������� ������ � �������� �� ����� ����������)

% INPUT DATA (� ��� �� ����� ��������� ��������� �������)
init_speed = boat_speed;%�/� (��������� ��������  �� ������)

% �������� ���� � �����
k = boat_k;%        (����������� ��������� �������� ���� �����)
k1 = boat_k1;%       (����������� ������������� ����)
vb = boat_vb;%�/�     (�������� �����)

% �����
m = boat_m;%��    (����� �����)
S = boat_S;%�^2     (������� ������)

% ��������� ���������� ��� ������� ����
curr_time = time;%                         ������� �������� ������ � ������ ��������
curr_gamma = acot( (init_x-dest_x)/(init_y-dest_y) );% ������� ���� �������� �����. 
curr_pos = boat_pos;
curr_speed = init_speed;%                ������� ��������, �������� ������ ��������

curr_x = init_x;% ��� ������� ��������� �� �
curr_y = init_y;% ��� ������� ��������� �� �

% ���� �� ������ �� ����� ���������� - ����� �������
while( dest_y - curr_y > 0  )

    % ���� (���������, ��� ����� ��������)
    av = 90*(pi/180);
    gamma = curr_gamma;    % (���� ����������� �������� ����� ������������ �� � ����)
    
    % ���� (�����������, ��� ����� ��������)
    theta = abs(av-gamma)/2;
    alpha = abs( abs(av-gamma) - theta ); % (���� ����� ������� � ������)
    
    
    % ���� curr_alpha > 90
    if (gamma < 0)
        gamma = pi + gamma;
        theta = abs(gamma - av)/2;
        alpha = theta;
    end
    
    
    % ���� (� ������ ��������)
    %delta =  acos((2*vb^2 - curr_speed^2)/(2*vb^2));    % (���� ���������� �������� ��� ������ ������������ ��������)
    %av =     90*(pi/180) - delta;    % (���� ����� ������������ ����)
    %gamma = curr_gamma;    % (���� ����������� �������� ����� ������������ �� � ����)
    %theta = (av-gamma)/2 - delta/2;         % (���� ����� ������� � ����� �����)
    %if(theta < 0)
    %    theta = 0;
    %end
    %alpha = theta;
    
    % ��������� ������� 
    x10 = curr_pos;%�       (��������� ��������� �� ������)
    x20 = curr_speed;%�/�     (��������� ��������  �� ������)
    tspan  = [curr_time curr_time + Tdelay];%� (���������� �������)
    x0 = [x10 x20];

    % ���������
    opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
    [t,z] = ode45(@(t,z) odefcn(t,z,k1,vb,alpha,theta,k,S,m), tspan, x0,opts);

    % ������� � ����� ������� ��������� (����)
    x0 = curr_x;
    y0 = curr_y;
    x = z(:,1)*cos(gamma)+x0;
    y = z(:,1)*sin(gamma)+y0;  

    % ��������� ����� ��������
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

