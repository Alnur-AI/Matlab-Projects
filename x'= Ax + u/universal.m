clc
clear

% ������ ��������� ��������� ������������ ����������� ����������
% ��� ��������� x'= Ax + u, x - ������ ������� �������
% �� ������� ������� [t0, T]
% M_0 - ��������� ���������
% U - ������� ����������
% ����������: J = min (t - t0)

%% ����� 1: ������� ����������� ��������, ������� A � ������� �������
% ����� ����� ������ ��������� ��������� F ����� �������� ����������.
% ������ ����� ���� ��������� F ����� � (0,0)
% ���� �� ��� �� ���, ����� ���� ����������� � ����� � ������ ����������.
% ��� �� �������� ��� 
% theta = linspace(0,2pi) - ���� �� 0 �� 2pi
% r = r(theta) - ������
% f0 = (x0,y0) - ���������� ������ ���� 
% ���������, ��� ����� (x,y), ������� �� ������� F ����� ���������
% ��������� �������: 
% x = r*cos(theta) + x0;
% y = r*sin(theta) + y0;
% ����������: ������� r(theta) ������ � ��������� �������� r_F,  

clc
clear

% �������
A = [1,1;-1,0];

% ��������� ����������
N_t = 70;%��������� �� ��������� t
t0 = 0;
T = 1;
t = linspace(t0,T,N_t);

% ������ psi
N_psi = 70;% ��������� �� ��������� psi, M_0, U, M_1
theta = linspace(0,2*pi,N_psi);

% ��������� ���������
cM_0 = [1;0];% ���������� ������ ���� ��������� M_0

i = 1:N_psi;
r = zeros(N_psi,1);
r(i) = r_M_0(theta(i));

M_0 = zeros(N_psi,2);% ����� �� ����� ��� ����� �� ������� M_0
M_0(:,1) = r(:).*cos(theta(:)) + cM_0(1);
M_0(:,2) = r(:).*sin(theta(:)) + cM_0(2);

% ��������� ����������
cU = [0;0];% ���������� ������ ���� ��������� U

i = 1:N_psi;
r = zeros(N_psi,1);
r(i) = r_U(theta(i));

U = zeros(N_psi,2);% ����� �� ����� ��� ����� �� ������� U
U(:,1) = r(:).*cos(theta(:)) + cU(1);
U(:,2) = r(:).*sin(theta(:)) + cU(2);

% ������ ��������� M_0 � U
figure (1)
hold on
plot(M_0(:,1),M_0(:,2),'-*');
plot(U(:,1),U(:,2),'-*');
legend('M_0','U')
hold off

clear i r

%% ����� 2: ���������� ������� ������� C(X(T),psi)
% ������ ����� ����� ��������� ������������ �� ��������� ������:
% CX(psi,t) = C(x(t),psi) = C(M_0, expm((t-t0)*A')*psi) + integral(   C(U,expm((t-s)*A')*psi )    ,t0,t,ds)

CX = zeros(N_psi,N_t);

% ��� ������� ���������� ������� psi(theta)
for j = 1:N_psi
    psi_j = [cos(theta(j)); sin(theta(j))];
    
    % ��� ������ ���������� t(i)
    for i = 1:N_t
        
        % ��������� �������� �� t0 �� t
        % integral(t0,t, C(U, expm((t-s)*A')*psi ),  ds );
        %%{
        CU = zeros(i,1);
        for k = 1:i
            CU(k) = C(U, expm((t(i)-t(k)).*transpose(A))*psi_j);
        end
        integral_CX = trapz(t(k),CU(k),2);
        %%}
        
        
        % ��������� �������� �� t0 �� T
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

%% ����� 3.1: ������ ������ CX � ���������� ����������� CX(theta,t)
% ��� ��� �� ����� ���� ����� ������� ������� C(x(t),psi)
% � ���������� �����, �� ����� ��������� ���������� ������
% �� ��� Ox - �����, Oy - ����, Oz - �������� ������� ������� � ���� �����
[XX,YY] = meshgrid(t,theta);

figure (1);
surf(XX,YY,CX)
xlabel('Time')
ylabel('psi(Angle)')
zlabel('CX')

clear XX YY

%% ����� 3.2: ������ ������ CX � �������������� ����������� CX(theta,t)
% �� ������ �������� ��� ������� CX(psi,t) ����� ������ 2pi �� ��� Ox
% �������������, ����������� ��� ������� � 
% �������������� ����������� ����� ���������
% ����� �������
% � ����� ������:
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

%% ����� 4.1: ������ x1(t) x2(t) �� �������:
%c = CX;% ����� c_i = CX(psi_i,x);
% x1(i)*cos(theta(j))   + x2(i)*sin(theta(j)) = c(j,i)
% x1(i)*cos(theta(j+1)) + x2(i)*sin(theta(j+1)) = c(j+1,i)
% ����� ������� �� ������������ ���������� (x1(t),x2(t))
% � ������������ �������� ������������ �� t0 �� T
% ������� ��������� ����������:
% xfull(psi,t) � ����� ������ ��� ������������ ����� 
% ����� ������������ �� t0 �� T


xfull = zeros(N_psi,N_t,2);

for j = 1:N_psi-1
   for i = 1:N_t
       xfull(j,i,:) = [cos(theta(j)) sin(theta(j)) ; cos(theta(j+1)) sin(theta(j+1))  ]\[CX(j,i);CX(j+1,i)];
   end
end
i = 1:N_t;
xfull(N_psi,i,:) = xfull(1,i,:);

%% ����� 4.2: �������� M_1 - ��������� ������������
M_1 = zeros(N_psi,2);

M_1(:,1) = xfull(:,end,1);
M_1(:,2) = xfull(:,end,2);

%% ����� 5.1: ������ � 2D ������������ ���������� (x1(t),x2(t)) �� M_0 �� M_1   
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

%% ����� 5.2: ������ � 3D ������������ ���������� (x1(t),x2(t),t) �� M_0 �� M_1  

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

%% ����� 6.1: ������ � 2D M_0 � M_1
figure (7)
hold on
    plot( M_0(:,1), M_0(:,2), '-*', mean(M_0(:,1)), mean(M_0(:,2)), '*' );
    plot( M_1(:,1), M_1(:,2), '-*', mean(M_1(:,1)), mean(M_1(:,2)), '*' );
    legend('M_0','����� ���� M_0','M_1','����� ���� M_1');
hold off





