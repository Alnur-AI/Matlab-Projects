%% Static points
b = 1.422;
c = -0.109;
f = @(n,a) n.*exp((b./n)-c./(n.^2)-a)
k = 100;
x_1 = -3;
x_2 = 3;
set_a = [x_1:1/k:x_2];
q = x_1:1/k:x_2;
plot(set_a,f(set_a,1),q,q);
axis([x_1 x_2 x_1 x_2])
legend({'f(N)','Diag'}, 'Location', 'northwest');
xlabel('a') 
ylabel('f(N)')
%% N_1
b = 1.422;
c = -0.109;
k = 10;
sys_1 = @(a) ((b^2 + sqrt(b^2-4*c.*a)-2*c.*a)./(2.*a));
sys_2 = @(a) ((b + sqrt(b^2-4*c.*a))./(2.*a) - 2*c/b);
axis_s = @(a) 0.*a;
linsp_a = [-4.638:1/k:5];
gcf = plot(linsp_a, sys_1(linsp_a), linsp_a, sys_2(linsp_a), linsp_a, axis_s(linsp_a));
legend({'Ineq 1','Ineq 2', '0-axis'}, 'Location', 'northwest');
xlabel('a') 
ylabel('Ineq 1, Ineq 2')

%% N_2
b = 1.422;
c = -0.109;
k = 10;
sys_1 = @(a) ((b^2 - sqrt(b^2-4*c.*a)-2*c.*a)./(2.*a));
sys_2 = @(a) ((b - sqrt(b^2-4*c.*a))./(2.*a) - 2*c/b);
axis_s = @(a) 0.*a;
linsp_a = [-4.638:1/k:5];
plot(linsp_a, sys_1(linsp_a), linsp_a, sys_2(linsp_a), linsp_a, axis_s(linsp_a));
legend({'Ineq 1','Ineq 2', '0-axis'}, 'Location', 'northwest');
xlabel('a') 
ylabel('Ineq 1, Ineq 2')

%% Bif Diagram
b = 1.422;
c = -0.109;
f = @(n,a) n.*exp((b./n)-c./(n.^2)-a);
n_0 = 1;
k_1 = 2700;
k_2 = 1500;
a_0 = 0.2;
delta = 0.001;
plot(0,0);
hold on
for i = 1:1:k_1
    x = zeros(k_2);
    t = n_0;
    for j = 1:1:k_2
        q = f(t,a_0);
        t = q;
        x(j) = t;
    end
    y = x(k_2-200:k_2);
    A=repmat(a_0,201,1);
    plot(A,y,'ko','MarkerSize',1);
   % axis([2.73 2.8 0 40])
    a_0=a_0+delta;
end
xlabel('a') 
ylabel('N_t') 


%% Cycles
clc
clear
a=2.73:0.001:2.8;
b = 1.422;
c = -0.109;
f = @(x,y)x.*exp(b./x-c./(x.^2)-y)
f1 = @(x, y) f(x, y).*exp(b./f(x, y)-c./(f(x, y).^2)-y);
f2i = @(x) f1(x, 2.763) - x;
t = fzero(f2i, 21);
%f2 = @(x, y) f1(x, y).*exp(b./f1(x, y)-c./(f1(x, y).^2)-y);
%f3 = @(x) f2(x, 2.763) - x;
%t = fzero(f3, 21);
disp(t);
disp(f(t, 2.763));
disp(f(f(t, 2.763), 2.763));

%% Cycles graphs
clc
clear
b = 1.422;
c = -0.109;
f = @(x)x.*exp(b./x-c./(x.^2)-2.763)
x = 0.08:0.01:23;
f1 = f(x);
figure();
axis([0 23 0 23]);

plot(x, f1,'r', x, x, 'b');
xlabel('N');
ylabel('f(N), N');
hold on;
axis([0 23 0 23]);
u = 2.763;
p0 = [u f(u)];
p1 = [f(u) f(u)];
line(p0, p1);
hold on;
p0 = [f(u) f(u)];
p1 = [f(u) u];
line(p0, p1);
hold on;
p0 = [f(u) u];
p1 = [u u];
line(p0, p1);
hold on;
p0 = [u u];
p1 = [u f(u)];
line(p0, p1);
hold on;
plot(u, f(u), 'm*', f(u), u, 'g*');
hold on;
u = 21.19;
p0 = [u f(u)];
p1 = [f(u) f(u)];
line(p0, p1);
hold on;
p0 = [f(u) f(u)];
p1 = [f(u) u];
line(p0, p1);
hold on;
p0 = [f(u) u];
p1 = [u u];
line(p0, p1);
hold on;
p0 = [u u];
p1 = [u f(u)];
line(p0, p1);
hold on;
plot(u, f(u), 'm*', f(u), u, 'g*');
hold off;

%% Lyapunov
clc
clear
b = 1.422;
c = -0.109;
f = @(x, y) log(abs((1-b./x + 2*c./(x.^2)).*exp(b./x-c./(x.^2)-y)));
axis_s = @(a) 0.*a;
n = 1000;
dx = 0.01;
linspx = [0:dx:3];

a = -1;

plot(linspx, f(linspx,a), linspx, axis_s(linspx));
axis([0 3 -7 7]);
xlabel('a') 
ylabel('Lyapunov Index') 



