%%Задание 1
a = 0;
b = 1;
n = 5;
space1 = linspace (a,b,n);
space2 = linspace (a,b,2*n);
g = @(x) exp(-(x.^2));
compareInterp(space1,space2,g);
%% Задание 2
%% Пример 1
% "nearest"
%  + Самые быстрые вычисления и минимальный расход памяти (ТОП1)
%  - функция получается разрывной
%  - довольная неточная интерполяция при растущей функции
a = 5;
b = 10;
n = 10;
space1 = linspace (a,b,n);
space2 = linspace (a,b,3*n);
g = @(x) x.^2;
v = g(space1);
vn = interp1(space1,v,space2,"nearest");
vs = interp1(space1,v,space2,"spline");
vv = g(space2);
plot(space2,vv,space2,vn)
disp('Максимальная погрешность метода ближ. соседа');
disp(abs(max(vv-vn)));
disp('Максимальная погрешность метода сплайнов ');
disp(abs(max(vv-vs)));
title("Недостатки метода ближайшего соседа");
legend("Функция","Интерполяция");
%% Пример 2
% "linear"
% + Обладет неприрывностью
% ++- Относительно быстрые вычисления и расход памяти (ТОП2)
% - если у функции большая вторая производная, то погрешность может быть неприемлимая
a = 2*pi;
b = 10*pi;
n = 9;
space1 = linspace (a,b,n);
space2 = linspace (a,b,5*n);
g = @(x) 10*sin(x./2).*sin(x./3);
v = g(space1);
vl = interp1(space1,v,space2,"linear");
vs = interp1(space1,v,space2,"spline");
vv = g(space2);
plot(space2,vv,space2,vl,space2,vs)
disp('погрешность линейного метода');
disp((trapz(space2,(vv-vl).^2)).^(1/2));
disp('погрешность метода сплайнов ');
disp((trapz(space2,(vv-vs).^2)).^(1/2));
title("Недостатки метода линейного");
legend("Функция","Интерполяция линейная","Интерполяция сплайнами");
%% Пример 3
%" cubic"
% + Класс неприрывности С1
% --+ довольно сложные вычисления (ТОП3)
% -для построения нужно 4 точки 
a = 0;
b = 3;
n = 1000;
space1 = linspace (a,b,n);
space2 = linspace (a,b,3*n);
g = @(x) sin(x.^3)+x.^3;
v = g(space1);
vl = interp1(space1,v,space2,"linear");
vc = interp1(space1,v,space2,"cubic");
vv = g(space2);
plot(space2,vv,space2,vl,space2,vc)
disp('Максимальная погрешность метода линейного');
disp(abs(max(vv-vl)));
disp('Максимальная погрешность метода кубов ');
disp(abs(max(vv-vc)));
title("Недостатки метода ближайшего соседа");
legend("Функция","Интерполяция соседом","Интерполяция кубами");
%% Пример 4 
% "spline"
a = 5;
b = 10;
n = 5000000;
space1 = linspace (a,b,n);
space2 = linspace (a,b,3*n);
g = @(x) 0.000001*x;
v = g(space1);
t1 = 0;
t2 = 0;
z = 10;
for i = 1:z
    tic();
    vn = interp1(space1,v,space2,"nearest");
    t1 = t1+toc();
    tic();
    vs = interp1(space1,v,space2,"spline");
    t2 = t2+toc();
end
disp('Время метода ближ соседа')
disp(t1/z);
disp('Время метода сплайнов')
disp(t2/z);
vv = g(space2);
disp('Максимальная погрешность метода ближ. соседа');
disp(abs(max(vv-vn)));
disp('Максимальная погрешность метода сплайнов ');
disp(abs(max(vv-vs)));
% + Класс неприрывности С2
% --- довольно сложные вычисления (ТОП4)
% -для построения нужно 4 точки 
%% Пример 5
a = 0;
b = 1;
n = 3;
space1 = linspace (a,b,n);
space2 = linspace (a,b,3*n);
g = @(x) x.^2;
v = g(space1);
vn = interp1(space1,v,space2,"lirear");
vv = g(space2);
plot(space2,vv,space2,vn)
title("Бонус интрерполяции линеным методом");
legend("Функция","Интерполяция");
%% Задание 3
h_st = 0.0001;
h_en = 0.01;
n = 1000;
h = linspace(h_st,h_en,n);
delta_g = h.*0;
delta_f =h.*0;
g = @(x) x.^4;
f = @(x) sin(20*x);
for i = 1:length(h)
    space1 = linspace (1,2,round(1/h(i)));
    space2 = linspace (1,2,round(2/h(i)));
    v_g = g(space1);
    v_f = f(space1);
    vs_g = interp1(space1,v_g,space2,"spline");
    vs_f = interp1(space1,v_f,space2,"spline");
    vv_g = g(space2);
    vv_f = f(space2);
    delta_g(i) = abs(max(vv_g-vs_g));
    delta_f(i) = abs(max(vv_f-vs_f));
end
maxp_g = 24;
maxp_f = 160000;
aprior_g = 5/384*maxp_g*h.^4;
aprior_f = 5/384*maxp_f*h.^4;

%plot(h,delta_g,h,aprior_g);
%legend("post x.^4","aprior x.^4");
plot(h,delta_f,h,aprior_f);
legend("post sin(10*x)","aprior sin(10*x)");

title('Interpol');
xlabel('h - шаг сетки');
%% Задание 4
a = 0;
b = 1;
k = 100;
% поточечная, но не равномерная 
%fn=@(x,n) x.^n;
%f = @(x) (x>=0)*0+(x==1)*1;
%convergenceFunc(fn,f,a,b,k,'Pointwise');
% сходимость в среднем, но не поточечная и не равномерная
fn  = @(x,n) 1*((x<= (1+mod(n,2.^(floor(log2(n)))))*2.^(-floor(log2(n)))) & (x>= mod(n,2.^(floor(log2(n))))*2.^(-floor(log2(n)))));
f = @(x) 0*(x<=1);
convergenceFunc(fn,f,a,b,k,'Average convergence');
%% Задание 5
f = @(x) x+sin(x)+cos(x);
a = -3;
b = 3;
n = 50;
%fourierApprox(f,a,b, n, "fourier");
%fourierApprox(f,a,b, 12, "hermite");
%fourierApprox(f,a,b, 10, "laguerre");
%% Задание 7
f = @(t) 2*sin(2*t+2);
g = @(t) 1*sin(t);
N = 5;
t0 = 0; 
t1 = pi-0.1;
 
mas = getEqual(f,g,t0,t1,N);

%% Задание 8

%% Задание 9

%% Задание 10



%% Задание 14
s = struct('ax',-10,'bx',10,'ay',-10,'by',10,'az',-10,'bz',10,'n',100,'fc','red','ec','None','radius',20,'light',1,'transparence',1);
drawBall(2,s);
%% Задание 15
    alphas = [0.5,1,2,1,2,1,2,1,2,1,2];
    colors=['b' 'r' 'y' 'r' 'y' 'r' 'y' 'r' 'y' 'r' 'y'];
    edges = ['b' 'r' 'y' 'r' 'y' 'r' 'y' 'r' 'y' 'r' 'y'];
    drawManyBalls(alphas,colors,edges);
%% Задание 16
function res = convergenceFunc(fn,f,a,b, n, convType)
mov(1:n) = struct('cdata', [],'colormap', []);
x = linspace(a,b,300);
y = f(x);
for i = 1:n
    yn  = fn(x,i); 
    met1 = abs(max(yn-y));
    met2 = abs(trapz(x,yn-y));
    plot(x,yn,x,y); % Не используется hold on: каждый новый кадр затирает предыдущий
    str = string(convType)+' Сходимость абсолютная '+string(met1)+' Сходимость в среднем ' + string(met2);
    title(str);
    legend("fn","function");
    mov(i) = getframe();
end
movie(mov);
end
function res = compareInterp(x,xx,f)
    v = f(x);
    vn = interp1(x,v,xx,"nearest");
    vl = interp1(x,v,xx,"linear");
    vs = interp1(x,v,xx,"spline");
    vc = interp1(x,v,xx,"cubic");
    plot(x,v,xx,vn,xx,vl,xx,vs,xx,vc);
    %xx,vn,xx,vl,xx,vs,xx,vc
    legend("function","nearest", "linear", "spline", "cubic");
    title('Interpol');
end
function res = drawBall(alpha,params)
params.radius = params.radius.^2;
[X,Y,Z] = meshgrid(params.ax:(params.bx - params.ax)/params.n:params.bx,params.ay:(params.by - params.ay)/params.n:params.by,params.az:(params.bz - params.az)/params.n:params.bz);
if (isinf(alpha))
    f = @(x,y,z) max(max(abs(x),abs(y)),abs(z));
else
    f = @(x,y,z) abs(x).^(alpha)+abs(y).^(alpha)+abs(z).^(alpha);
end
F = f(X,Y,Z);
g = isosurface(X,Y,Z,F,params.radius);
disp(g);
if (~isempty(g.faces))
    p = patch(g);
    p.FaceColor = params.fc;
    p.EdgeColor = params.ec;
    p.FaceAlpha = params.transparence;
    if (params.light ~= 0)
        camlight;
    end
    view(3);
else
    disp('Множество пустое');
end

end
function res = drawManyBalls(alphas,colors,edges)
    transp = 0.1:(0.9/length(alphas)):1;
    s = struct('ax',-1,'bx',1,'ay',-1,'by',1,'az',-1,'bz',1,'n',10,'fc',colors(1),'ec',edges(1),'radius',1,'light',1,'transparence',1-transp(1));
    drawBall(alphas(1),s)
    hold on
    for i = 2:length(alphas) 
        s = struct('ax',-1,'bx',1,'ay',-1,'by',1,'az',-1,'bz',1,'n',10,'fc',colors(i),'ec',edges(i),'radius',1,'light',0,'transparence',1-transp(i));
        drawBall(alphas(i),s)
    end
    hold off
end
function fourierApprox(f,a,b, n, meth)
    if (string(meth) == "fourier")
        fourierApp(f,a,b, n);
    end
    if (string(meth) == "hermite" )
        hermiteApp(f,a,b, n);
    end
    if (string(meth) == "laguerre")
        laguerreApp(f,a,b, n);
    end
end
function fourierApp(f,a,b, n)
    
    x = linspace(a,b,2000);
    x_s = x;
    y = f(x);
    x = pi*(x-a)/(b-a);
    sumfourier = zeros(n+1,length(x));
    for i = 1:n
        g_n = Getfun1(i);
        y_n = g_n(x);
        a_n = trapz(x,y_n.*y);
        sumfourier(i+1,:)=sumfourier(i,:)+2*a_n.*y_n/pi;
    end
    
    max_four  = max(max(sumfourier));
    min_four = min(min(sumfourier));
    mov(1:n+1) = struct('cdata', [],'colormap', []);
    for i = 1:n+1
        plot(x_s,y,x_s,sumfourier(i,:));
        axis([a b min_four max_four])
        legend("function","fourierApprox");
        title('fourierApprox');
        mov(i) = getframe();
    end
    movie(mov);
end
function res = Getfun1(n)
    res = @(x) sin(n*x);
end
function hermiteApp(f,a,b, n)
    x = linspace(a,b,100);
    y = f(x);
    sumfourier = zeros(n+1,length(x));
    for i = 0:n-1
        g = Getfun2(i);
        g_s = g(x);
        x_n = linspace(-100,100,2000);
        g_n = g(x_n);
        y_n = f(x_n);
        exp_n = exp(-(x_n).^2);
        a_n = trapz(x_n,g_n.*y_n.*exp_n);
        norm_n = trapz(x_n,g_n.*g_n.*exp_n);
        a_n = a_n/norm_n;
        sumfourier(i+2,:)=sumfourier(i+1,:)+a_n.*g_s;
    end
    %рисум
    max_four  = max(max(sumfourier));
    min_four = min(min(sumfourier));
    mov(1:n+1) = struct('cdata', [],'colormap', []);
    for i = 1:n+1
        plot(x,y,x,sumfourier(i,:));
        axis([a b min_four max_four])
        legend("function","hermiteApprox");
        title('hermiteApprox');
        mov(i) = getframe();
    end
    movie(mov);
    
end
function res = Getfun2(n)
    res = @(x) hermiteH(n,x);
end
function laguerreApp(f,a,b, n)
    x_s = linspace(a,b,100);
    x = x_s-a;   
    y = f(x);
    sumfourier = zeros(n+1,length(x));
    for i = 0:n-1
        g = Getfun3(i);
        g_s = g(x);
        x_n = linspace(0,100,2000);
        g_n = g(x_n);
        y_n = f(x_n);
        exp_n = exp(-(x_n));
        a_n = trapz(x_n,g_n.*y_n.*exp_n);
        norm_n = trapz(x_n,g_n.*g_n.*exp_n);
        a_n = a_n/norm_n;
        sumfourier(i+2,:)=sumfourier(i+1,:)+a_n.*g_s;
        
    end
    %рисум
    max_four  = max(max(sumfourier));
    min_four = min(min(sumfourier));
    mov(1:n+1) = struct('cdata', [],'colormap', []);
    for i = 1:n+1
        plot(x_s,y,x_s,sumfourier(i,:));
        axis([a b min_four max_four])
        legend("function","LaguerreApprox");
        title('LaguerreApprox');
        mov(i) = getframe();
    end
    movie(mov);    
end
function res = Getfun3(n)
    res = @(x) laguerreL(n, x);
end
function res = getEqual(f,g,t0,t1,N)
n = 3000;
sp = linspace (t0,t1,n);
x = f(sp);
y = g(sp);
m = cat (2,x',y');
D = pdist(m);
Z = squareform(D);
dist = max(max(Z))/2;

ans_mas = 1:N-1;
ans_mas(end+1) = n;
dist_mas = zeros(1,N-1);

for zz = 1:(length(ans_mas)-1)
    dist_mas(zz)=abs(((f(sp(ans_mas(zz)))-f(sp(ans_mas(zz+1)))).^2+(g(sp(ans_mas(zz)))-g(sp(ans_mas(zz+1)))).^2).^(1/2));
end
while(true)

    [d_max,di_max] = max(dist_mas);
    if (di_max==1)
       disp('a');
       break;
    else
       ans_mas(di_max)=ans_mas(di_max)+1;
    end
    dist_mas(di_max)=abs(((f(sp(ans_mas(di_max)))-f(sp(ans_mas(di_max+1)))).^2+(g(sp(ans_mas(di_max)))-g(sp(ans_mas(di_max+1)))).^2).^(1/2));
    
    if (di_max>=2)
        dist_mas(di_max-1)=abs(((f(sp(ans_mas(di_max-1)))-f(sp(ans_mas(di_max)))).^2+(g(sp(ans_mas(di_max-1)))-g(sp(ans_mas(di_max)))).^2).^(1/2));

    end
    if ((max(dist_mas)-min(dist_mas))<0.05)
        break;
    end
end
disp('Distanse');
for zz = 1:(length(ans_mas)-1)
    disp(((f(sp(ans_mas(zz)))-f(sp(ans_mas(zz+1)))).^2+(g(sp(ans_mas(zz)))-g(sp(ans_mas(zz+1)))).^2).^(1/2))
end
disp('Distanse rr');
rr = linspace(t0,t1,N);
for zz = 1:N-1
    disp(((f(rr(zz))-f(rr(zz+1))).^2+(g(rr(zz))-g(rr(zz+1))).^2).^(1/2));
end

plot(x,y,f(sp(ans_mas)),g(sp(ans_mas)),'r*',f(rr),g(rr),'b*');
hold on
for k=1:length(ans_mas)
    text(f(sp(ans_mas(k))), g(sp(ans_mas(k))), num2str(k))
end
plot (f(t0),g(t0),'bp',f(t1),g(t1),'ko','MarkerSize',15);
hold off

res = ans_mas;
end