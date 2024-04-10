function [curr_x, curr_y, curr_pos, curr_speed, curr_time] = transfer_between_dots(init_x,init_y,dest_x,dest_y,boat_k,boat_k1,boat_vb,boat_m,boat_S, boat_pos,  boat_speed,time)

% ��������� ��������
Tdelay = 1;%� (��������� ���������� ��� �������� ������ � �������� �� ����� ����������)

% INPUT DATA (� ��� �� ����� ��������� ��������� �������)
init_speed = boat_speed;%�/� (��������� ��������  �� ������)

% �������� ���� � �����
k = boat_k;%        (����������� ��������� �������� ���� �����)
k1 = boat_k1;%       (����������� ������������� ����)
vb = boat_vb;%�/�     (�������� �����)
av = pi/2;

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

    
    % ���� (� ������ ��������)
    angle_CBA = pi - (av - curr_gamma);
    
    % � ��������������
    vb = sqrt(boat_vb^2 + curr_speed^2 -2*boat_vb*curr_speed*cos(angle_CBA) );
    delta = acos( (boat_vb^2 + vb^2-curr_speed^2)/(2*boat_vb*vb) );
    
    % ��� ������������
    %vb = boat_vb;
    %delta = 0;
    
    av =     90*(pi/180) + delta ;     % (���� ����� ������������ ����)
    gamma = curr_gamma;    % (���� ����������� �������� ����� ������������ �� � ����)
    theta = (av-gamma)/2 -delta/2;         % (���� ����� ������� � ����� �����)
    alpha =  theta;
    fprintf('vb = %f, av = %f, theta = %f\n',vb,av*180/pi,theta*180/pi);
    %disp((pi/2 - 2*theta - gamma+delta)*180/pi);
    
    % ���� curr_alpha > 90
    if (gamma < 0)
        av =     90*(pi/180) - delta;
        gamma = pi + gamma;
        theta = abs(gamma - av)/2;
        alpha = theta;
        disp((gamma)*180/pi);
    end
    
    
    
    %{
    % ���� (���������, ��� ����� ��������)
    av = 90*(pi/180);
    gamma = curr_gamma;    % (���� ����������� �������� ����� ������������ �� � ����)
    theta = abs(av-gamma)/2;
    alpha = abs( abs(av-gamma) - theta ); % (���� ����� ������� � ������)
    
    
    % ���� curr_alpha > 90
    if (gamma < 0)
        gamma = pi + gamma;
        theta = abs(gamma - av)/2;
        alpha = theta;
    end
    %}
    
    
    
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

