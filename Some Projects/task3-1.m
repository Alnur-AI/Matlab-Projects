%% Задание_1
% e-1 = 1/1! + 1/2! + ... 1/n! + a/n*n! где 0<a<1
% e-1 = 1/1! + 1/2! + ... 1/n! + 1/(n+1)! ...
% a/n*n! = 1/(n+1)! + 1/(n+2)! ...
% Sn = 1/(1*1!)+1/(2*2!)+1/(3*3!) + ... + 1/(n*n!)+ r(n)
% r(n) = 1/((n+1)*(n+1)!)+1/((n+2)*(n+2)!)+ ... 
% r(n)<= a/n*n! <= 1/(n*n!)
N = 15;
vec1 = ones(1,N);% [1 1 1 1 1 ...]
vec2 = cumsum(vec1);% [1 2 3 4 5 ... n]
vec3 = factorial(vec2);% [1 2 6 34 120 720 ... n! ]
pog = (vec1./vec2)./vec3; %[1 1/2*2! ...1/n*(n)! ]
y = cumsum(pog); % [S1 S2 S3 .... Sn]
y = abs(y(N)- y);
x = linspace(1,N-1,N);
plot(x,y,x,pog);
xlabel('N');
legend('abs(S-S_N)','1/(n*n!)');
%%  Задание_2
f =@(x) x.^(1/2);
g =@(x) tan(x);
F =@(x) x.^(1/2) - tan(x);
N = 100;
X = linspace(0,10*pi,N);
y_f = f(X);
y_g = g(X);
plot(X,y_f,X,y_g);
m = 7;
[cord_x,cord_y] = ginput(m);
hold on
for i = 1:m
    g = round((cord_x(i))/pi);
    
    solution = fzero(F,g*pi);
    plot (solution,f(solution),'r.','MarkerSize',30);
    plot (cord_x(i),cord_y(i),'b.','MarkerSize',30);
    disp(i);
    disp(solution)
    legend('x.^(1/2) - tan(x)');
    disp(abs(F(solution)));
end
%% Задание_3
eps = 0.00001;
f = @(x) ((abs(x)).^(1/2)).*sin(x.^(-2));
a = 10;
n = 100;
x = linspace(-a,a,n);
x_s_m = zeros(1,n);
for i = 1:n
    if (abs(x(i))>=1)
        x_s_m(i) = fzero(f,abs(x(i))/x(i));
    else
    	x_s_m(i) = fzero(f,x(i));
    end
end
plot(x,x_s_m);
xlabel('x');
legend('Ближайший корень fzero');
%% Задание_4
alf = 0.2;
tspan = [0 10];
y0 = [1 2];
opts = odeset('Events',@ bounceEvents);
disp(opts);
[t,y,te,ye,ie] = ode45(@(t,y) odefcn(t,y,alf),tspan,y0,opts);
    if (size(te)~=0)
        tspan(1) = te;
    else
        tspan(1) = tspan(2);
    end
    
    if (size(ye)~=0)
    	y0(1) = 0.001*ye(2)/abs(ye(2));
    	y0(2)= -ye(2);
    end
plot(t,y(:,1),'b',t,y(:,2),'r');
xlabel('Time (t)');
hold on;
while(tspan(1)+0.01<tspan(2))
    [t,y,te,ye,ie] = ode45(@(t,y) odefcn(t,y,alf),tspan,y0,opts);
    disp(te);
    plot(t,y(:,1),'b',t,y(:,2),'r');
    if (size(te)~=0)
        tspan(1) = te;
    else
        tspan(1) = tspan(2);
    end
    if (size(ye)~=0)
    	y0(1) = 0;
    	y0(2)= -ye(2);
    end
end
legend('x(t)','v(t)');
hold off 

%% Задание_5
tspan = [0 5];
y0 = [0 1 0 1 4 0 1 1 0 3 0 3];
%y0 = [0 1 0 1 4 0 1 1 0 3 0 0];

%y0 = [0 1 0 1 1 0 1 1 0 1 0 0];

m1 = 2;
m2 = 1; 
[t,y] = ode45(@(t,y) odef_5(t,y,m1,m2),tspan,y0);
m = length(t);

sum_x2 = sum(y(:,1).^2)+sum(y(:,7).^2);
sum_y2 = sum(y(:,2).^2)+sum(y(:,8).^2);
sum_z2 = sum(y(:,3).^2)+sum(y(:,9).^2);
sum_x = sum(y(:,1))+sum(y(:,7));
sum_y = sum(y(:,2))+sum(y(:,8));
sum_z = sum(y(:,3))+sum(y(:,9));
matrix = [sum_x2 sum_y sum_z; sum_x sum_y2 sum_z; sum_x sum_y sum_z2 ];
P = pinv(matrix);
mat_left = [-1; -1; -1];
coff = P*mat_left;
disp(coff);

max_X_1 = max(y(:,1));
max_X_2 = max(y(:,7));
max_X = max(max_X_1,max_X_2);

max_Y_1 = max(y(:,2));
max_Y_2 = max(y(:,8));
max_Y = max(max_Y_1,max_Y_2);

max_Z_1 = max(y(:,3));
max_Z_2 = max(y(:,9));
max_Z = max(max_Z_1,max_Z_2);

min_X_1 = min(y(:,1));
min_X_2 = min(y(:,7));
min_X = min(min_X_1,min_X_2);

min_Y_1 = min(y(:,2));
min_Y_2 = min(y(:,8));
min_Y = min(min_Y_1,min_Y_2);

min_Z_1 = min(y(:,3));
min_Z_2 = min(y(:,9));
min_Z = min(min_Z_1,min_Z_2);

x_vec = linspace(min_X,max_X,10);
y_vec = linspace(min_Y,max_Y,10);
[ set_x, set_y ] = meshgrid(x_vec,y_vec);
Z = (coff(1).*set_x+coff(2).*set_y)./coff(3);
if (coff(3) ==0)
    Z=0*set_x;
end
%coff(3) ==0
mov(1:m) = struct('cdata', [],'colormap', []);
for i = 1:m
    hold on
    surf(set_x,set_y,Z,'FaceAlpha',0.5);
    plot3(y(i,1),y(i,2),y(i,3),'b.','MarkerSize',30);
    plot3(y(i,7),y(i,8),y(i,9),'r.','MarkerSize',30);
    axis([min_X-1 max_X+1 min_Y-1 max_Y+1 min_Z-1 max_Z+1]);
    view(3);
    mov(i) = getframe();
    hold off
end

%% Задание_5_2
movie(mov);
%% Задание_6
tspan = [0 1];
xlabel('x');
ylabel('y')
hold on;
for i = -0.5:0.2:0.5
    for j = -0.5:0.2:0.5
        y0 = [i j];
        [t,y] = ode45(@(t,y) odef_1(t,y),tspan,y0);
        m =length(y(:,1));
        n = 7;
        m = fix(m/n)*n;
        norm = 0.2;
       for l = 1:n:m-n
                leng = (1/norm)*((y(l-1+n,1)-y(l,1)).^2+(y(l-1+n,2)-y(l,2)).^2).^(1/2);
                q = quiver(y(l,1),y(l-1+n,2),(y(l-1+n,1)-y(l,1))/leng,(y(l-1+n,2)-y(l,2))/leng,'Color','r');
                q.LineWidth = 3;
        end
        %q.AutoScaleFactor = 0.5;
        %q.AutoScale = 'off';
        %q.LineWidth = 0.8;
        %plot(y(:,1),y(:,2),'r');
    end
end
hold off

%% Задание_7
f1 = @(x,y) (x+y).^2+y.^4/2; 
f2 = @(x,y) x.*y;

tspan = [0 4];
hold on;
for i = -1:0.40:1
    for j = -1:0.40:1
        y0 = [i j];
        [t,y] = ode45(@(t,y) odef7_1(t,y),tspan,y0); 
        m =length(y(:,1));
        n = 5;
        norm = 0.1;
        m = fix(m/n)*n;
        for l = 1:n:m-n
            % f < 3 1*f1 
            % f2 < 1 3*f2
            if (f1(abs(y(l-1+n,1)),abs(y(l-1+n,2)))<1)
                leng = (1/norm)*((y(l-1+n,1)-y(l,1)).^2+(y(l-1+n,2)-y(l,2)).^2).^(1/2);
                q = quiver(y(l,1),y(l-1+n,2),(y(l-1+n,1)-y(l,1))/leng,(y(l-1+n,2)-y(l,2))/leng,'Color',[0.20+f1(abs(y(l-1+n,1)),abs(y(l-1+n,2)))/4 0 0]);
                q.LineWidth =1.5;
            end
        end
        %q = quiver(y(1:n:m,1),y(n:n:m,2),y(n:n:m,1)-y(1:n:m,1),y(n:n:m,2)-y(1:n:m,2),'r');
        %q.AutoScaleFactor = 0.6;
        %q.AutoScale = 'off';
        %plot(y(:,1),y(:,2),'r');
        %TODO цвета и длинну векторов
    end
end
X_ = linspace(-1,1,50);
Y_ = linspace(-1,1,50);
%X_ = linspace(0,1,50);
%Y_ = linspace(0,1,50);
[X,Y] = meshgrid(X_,Y_);
Z = f1(X,Y);
contour(X,Y,Z,25);
hold off
%% Задание_8
xmesh = linspace(0,1,5);
solinit = bvpinit(xmesh, @guess);
sol = bvp4c(@bvpfcn, @bcfcn, solinit);
hold on
f = @(t) exp(t)-2;
xmesh_2 = linspace(0,1,20);
y = f(xmesh_2);
y_1 = f(xmesh);
plot(sol.x,sol.y(2,:),xmesh_2,y,'r');
C_metrix = max(abs(y_1-sol.y(2,:)));
L2_metrix = (trapz(sol.x,(y_1-sol.y(2,:)).^2)).^(1/2);
disp('C_metrix');
disp(C_metrix);
disp('L2_metrix');
disp(L2_metrix);
legend('bvp4c','real solve');
hold off
%% Задание_9
x0 = [10 10];
%f = @(x,y) -(x+1).^2+y.^2+x.*y+x.^4;
%grad = {@(x) -2*x(1)-2+x(2)+4*x(1).^3, @(x) 2*x(2)+x(1)};
%gessian = {@(x) -2+12*x(1).^2, @(x) 1; @(x) 1, @(x) 2};
f = @(x,y) x.^2+y.^2-4;
grad = {@(x) 2*x(1), @(x) 2*x(2)};
gessian = {@(x) 2 , @(x) 0; @(x) 0, @(x) 2};
res = my_fminbnd(x0,grad,gessian);
x = -5:0.1:5;
y = -5:0.1:5;
[X,Y] = meshgrid(x,y);
Z = f(X,Y);
hold on
contour(X,Y,Z,'ShowText','on');
v = zeros(20,1);
for i=1:20
    v(i) = f(res(1,i),res(2,i));
    plot(res(1,i),res(2,i),'b.','MarkerSize',30);
    text(res(1,i)+0.1, res(2,i), num2str(i),'Color','r','FontSize',12)
end
plot(res(1,20),res(2,20),'r.','MarkerSize',30);
disp(res);
disp('my_val');
disp (v(20));
contour(X,Y,Z,v,'ShowText','on');
disp('fminbnd');
fun = @(x) x.^4-3*x.^3-5*x.^2+3;
[x_,fval] = fminbnd(fun,1,2*pi)
grad = {@(x) 4*x.^3-9*x.^2-10*x};
gessian = {@(x) 12*x.^2-18*x-10};
x0 = 5;
res_2 =my_fminbnd(x0,grad,gessian);
disp('my_fminbnd');
disp(res_2(20));
disp(fun(res_2(20)));
disp('abs');
disp(abs(res_2(20) - x_));
%% Задание_10
%% Функции
function dydt = odefcn(t,y,alf)
dydt = zeros(2,1);
dydt(1) = -sin(y(2))-alf*y(1);
dydt(2) = y(1);
end

function [value,isterminal,direction] = bounceEvents(t,y)
value = y(2)+1*(y(1)>=0);     % Detect height = 0
isterminal = 1;   % Stop the integration
direction = 0;   % Negative direction only
end
%Функция для задания 6
function dydt = odef_1(t,y)
dydt = zeros(2,1);
dydt(1) = y(1)-y(2);
dydt(2) = y(1)+y(2);
end
%Фокус
%dydt(1) = y(1)-y(2);
%dydt(2) = y(1)+y(2);
%Центр
%dydt(1) = -y(2);
%dydt(2) = y(1);
%Седло
%dydt(1) = -y(1);
%dydt(2) = 2*y(2);
%Узел 
%dydt(1) = 3*y(1);
%dydt(2) = 2*y(2);
%Д. Узел 
%dydt(1) = 2*y(1);
%dydt(2) = 2*y(2);
function dydx = bvpfcn(x,y)
dydx = zeros(2,1);
dydx = [y(1)
       y(1)];
end
function res = bcfcn(ya,yb)
res = [ya(2)+1
       yb(1)-2-yb(2)];
end
function g = guess(x)
g = [exp(x)
     exp(x)];
end
function dydt = odef7_1(t,y)
dydt = zeros(2,1);
dydt(1) = 2*y(2)-y(1)-y(2).^3;%x
dydt(2) = y(1)-2*y(2);%y
end
function dydt = odef7_2(t,y)
dydt = zeros(2,1);
dydt(1) = y(1)*y(2)-y(1).^3+y(2).^3;%x
dydt(2) = y(1).^2-y(2).^3;%y
end
function dydt = odef_5(t,y,m1,m2)
G = 6,674/100000000000;
dydt = zeros(12,1);
dydt(1) = y(4);%x1
dydt(2) = y(5);%y1
dydt(3) = y(6);%z1
dydt(4) = G*m2*(y(7)-y(1))/((y(1)-y(7)).^2+(y(2)-y(8)).^2+(y(3)-y(9)).^2).^(3/2);%vx1
dydt(5) = G*m2*(y(8)-y(2))/((y(1)-y(7)).^2+(y(2)-y(8)).^2+(y(3)-y(9)).^2).^(3/2);%vy1
dydt(6) = G*m2*(y(9)-y(3))/((y(1)-y(7)).^2+(y(2)-y(8)).^2+(y(3)-y(9)).^2).^(3/2);%vz1
dydt(7) = y(10);%x2
dydt(8) = y(11);%y2
dydt(9) = y(12);%z2
dydt(10) = G*m1*(y(1)-y(7))/((y(1)-y(7)).^2+(y(2)-y(8)).^2+(y(3)-y(9)).^2).^(3/2);%vx2
dydt(11) = G*m1*(y(2)-y(8))/((y(1)-y(7)).^2+(y(2)-y(8)).^2+(y(3)-y(9)).^2).^(3/2);%vy2
dydt(12) = G*m1*(y(3)-y(9))/((y(1)-y(7)).^2+(y(2)-y(8)).^2+(y(3)-y(9)).^2).^(3/2);%vz2
end
function res = my_fminbnd(x0,grad,gessian)
    N = 20;
    res = zeros(length(x0),N);
    res(:,1) = x0;
    for i=2:N
        X = zeros(length(x0),length(x0));
        for j =1:length(x0)
            for k=1:length(x0)
                q = gessian{j,k};
                g = q(res(:,i-1)');
                X(j,k) = g;
            end
        end
        X_ob = inv(X);
        grad_0 = zeros(1,length(x0));
        for k=1:length(x0)
            q = grad{k};
            grad_0(k) = q(res(:,i-1)');
        end
        grad_0 = grad_0';
        res(:,i) = res(:,i-1)-X_ob*grad_0;
    end
end
