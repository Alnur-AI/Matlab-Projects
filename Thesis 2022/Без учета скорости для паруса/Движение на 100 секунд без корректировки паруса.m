
clc
clear

% ���� � �����
k = 10;      %    (����������� ��������� �������� ���� �����)
k1 = 10;     %    (����������� ������������� ����)
vb = 5;%�/�      (�������� �����)

% �����
m = 1000;%��    (����� �����)
S = 10;%�^2     (������� ������)

% ���� (���������)
av =     90*(pi/180);    % (���� ����� ������������ ����)
gamma = 60*(pi/180);    % (���� ����������� �������� ����� ������������ �� � ����)
delta =  0*(pi/180);    % (���� ���������� ��������)

% ���� (�����������)
theta = (av-gamma)/2 + delta;         % (���� ����� ������� � ����� �����)
alpha = abs( abs(av-gamma) - theta ); % (���� ����� ������� � ������)

% ��������� ������� 
x10 = 0;%�       (��������� ��������� �� ������)
x20 = 0;%�/�     (��������� ��������  �� ������)
tspan  = [0 1000];%� (���������� �������)
x0 = [x10 x20];

% ���������
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
[t,z] = ode45(@(t,z) odefcn(t,z,k1,vb,alpha,theta,k,S,m), tspan, x0,opts);

% ������� � ����� ������� ��������� (����)
x0 = 0;
y0 = 0;
x = z(:,1)*cos(gamma)+x0;
y = z(:,1)*sin(gamma)+y0;

% ������ �������
%hold on
%plot(x,y);

%boat_scale = y(end)/8;
%quiver(x(end),y(end),boat_scale*cos(gamma),boat_scale*sin(gamma),'black');
%quiver(x(end),y(end),boat_scale*cos(theta+gamma),boat_scale*sin(theta+gamma),'b');% theta [0,180]
%quiver(x(end),y(end),boat_scale*cos(av),boat_scale*sin(av),'magenta');

%axis([-1 1 -1 1])
%plot(t,z(:,1),'-o',t,z(:,2),'-.')
plot(t,z(:,2),'-.')







%{
% ������ ������������� �������
A = -k1/m;
B = -k*S*(sin(theta)^2)/m;
C = vb*sin(alpha)*sin(theta)*k*S/m;

if (B^2-4*A*C>0)
    X_A = (-B + sqrt(B^2-4*A*C) )/(2*A);
    X_B = (-B - sqrt(B^2-4*A*C) )/(2*A);

    % ������� �������������
    f = @(t) (X_B-X_A)/(1 - abs(X_B/X_A)*exp(A*(X_B-X_A)*t) ) + X_A
    g = @(t) (X_B-X_A)/(1 + abs(X_B/X_A)*exp(A*(X_B-X_A)*t) ) + X_A % ��� av = 270, theta = [0,180)
end
if (0 <= B^2-4*A*C  && B^2-4*A*C <= 0.000000001 )
    disp(B^2-4*A*C)
end
if (B^2-4*A*C<0)
    disp(B^2-4*A*C)
end



for i = 1:size(t)
    an(i,1) =  f(t(i));
    an(i,2) =  g(t(i));
end
%}

x(end)
y(end)
theta*(180/pi)
alpha*(180/pi)







function dydt = odefcn(t,z,k1,vb,alpha,theta,k,S,m)
  dydt = zeros(2,1);
  dydt(1) = z(2);
  dydt(2) = ((-k1)*(z(2)^2) + ( vb*sin(alpha)-z(2)*sin(theta) )*sin(theta)*k*S)/m;
end

