% Вариант 1
%% task 1 (complete)
clear

% Input data
a = 0;
b = 10;
n = 10000;

% Plot function
x = linspace(a,b,n);
y = sin( log(1+abs(x)) + x.^2);

% Find all maxima and minima 
mini = islocalmin(y);
maxi = islocalmax(y);

%Plot function and all extremums
plot(x,y,x(mini),y(mini),'*',x(maxi),y(maxi),'*');

%% task 2 (complete)
clear

% Is natural?
n = 4;
isnatural(n);

% First
I = 9:18:n;
disp(I);

% Second
a = zeros(n);

a(:,n) = 1:n;
a = repmat (a(:,n),1,n);

% Third
B = zeros(n,n+1);
B(:) = 1:1:n*(n+1);
B = B';


c = (B(:))';
c = c';

D = B( : , (end-1):end );
%% task 3 (complete)
clear

% 3.1 Create matrix
a = 1 + floor((316-1).*rand(7));

% 3.2 Max on diagonal
maxV = max(diag(a));

% 3.3 ??????????????????

% 3.4 Order (by bubble sort)

B = sortrows(a);
disp(B);

%% task 4 (complete)
clear

x = [1 2 3];
y = [4 5 6];
a = (x') * y;

disp(x);
disp(y);
disp(a);
%% task 5 (complete)
clear

n = 7;
if isprime(n) == 0
    disp('not prime');
    return
end

a = rand(n);
b = rand(n,1);
if det(a) == 0
    disp('det equals zero');
    return
end

x1 = linsolve(a,b);
x2 = a\b;
x3 = inv(a)*b;

norm_x1x2 = norm(x1-x2);
norm_x2x3 = norm(x2-x3);
norm_x3x1 = norm(x3-x1);

%% task 6 (complete)
clear

m = 5;
n = 9;

a = rand(n,1);
b = rand(m,1);

maxv = max( max(a) - min(b), max(b) - min(a) );
%% task 7 (complete)
clear

n = 3;
k = 3;

a = [1 2 3; 5 0 2; 2 0 0];

B = zeros(n,n);
B1 = zeros(n,n);
B2 = zeros(n,n);
B3 = zeros(n,n);

B2 = a*(a'); % B2(i,j) = { a(i,:)*a(j,:) } 
temp = diag(B2); % temp(i,1) = { a(i,:)*a(i,:) }
B1 = repmat(temp,1,n); % B1(i,:) = { a(i,:)*a(i,:) }
B3 = B1'; % B1(:,j) = { a(:,j)*a(:,j) }

B = sqrt(B1-2*B2+B3);
%% task 8 (complete)
clear 

n = 3;

v = 0:1:2^n-1;

v = de2bi(v,'left-msb');


disp(v);

%% task 9 (complete)
clear

NUM = 200;
CHECK = 100;
ml_time = zeros(NUM,1);
my_time = zeros(NUM,1);
x = 1:NUM;
m_norm = zeros(NUM,1);



for num = 1:NUM
    c_ml_time = zeros(CHECK,1);
    c_my_time = zeros(CHECK,1);
    
    for count = 1:CHECK
        m = num;
        k = 10;
        n = num;

        a = rand(m,k);
        B = rand(k,n);

        tic;
        C = a * B;
        c_ml_time(count,1) = toc;

        tic;
        CC = my_multiply(a,B);
        c_my_time(count,1) = toc;

        m_norm(num,1) = norm (C-CC);
        disp(num)
    end
    ml_time(num, 1) = sum(c_ml_time)/CHECK;
    my_time(num, 1) = sum(c_my_time)/CHECK;
end
clear num;

plot(x, ml_time, x, my_time);

%% task 10 (complete)
clear

X = [nan,0,2,4;
    nan, 0,6,3;
    1,   5,nan,nan;
    1,   5,2,nan];

meanX = nanOutMean3(X);
disp(meanX);

%% task 11 (complete) 
clear

n = 100000;
a = 10;
sig = sqrt(21);

F = normrnd(a*ones(1,n), sig);


n_1sig = sum (a-1*sig <= F & F <= a+1*sig)*100/n;
n_2sig = sum (a-2*sig <= F & F <= a+2*sig)*100/n;
n_3sig = sum (a-3*sig <= F & F <= a+3*sig)*100/n;

Fnew = a-3*sig <= F & F <= a+3*sig;
Fnew = Fnew.*F;

Fnew = nonzeros(Fnew);
%% task 12 (complete)
clear

a = -pi;
b = pi;
n = 10000;

% Plot |h - h/2|
for i = 1:n
    
    h(i) = (b-a)/(2*i);
    
    xx = a : h(i) : b;
    yy = dSi(xx);
    
    
    
    Ihr = my_rectangles(xx,yy);
    Iht = trapz(xx,yy);
    Ihs = my_simpson(xx,yy);
    
    h2 = (b-a)/(4*i);
    xx2 = a : h2 : b;
    yy2 = dSi(xx2);
    
    Ih2r = my_rectangles(xx2,yy2);
    Ih2t = trapz(xx2,yy2);
    Ih2s = my_simpson(xx2,yy2);
    
    Ir(i) = abs(Ihr - Ih2r);
    It(i) = abs(Iht - Ih2t);
    Is(i) = abs(Ihs - Ih2s);
    
end

clear xx yy xx2 yy2 Ihr Iht Ihs Ih2r Ih2t Ih2s h2


figure(1);
loglog(h, Ir, h, It, h, Is);
title('Time difference');
legend({'|Ir-Ir2|','|It - It2|', '|Is - Is2|'},'Location','southwest')
xlabel('h');
ylabel('|dy - dy_n|');



%plot Si(x) by all methods

h = (b-a)/n;
x = a:h:b; %n+1 elements


ys1 = zeros(1,n+1);
ys2 = zeros(1,n+1);
ys3 = zeros(1,n+1);

for i = 1:n
    xx = a: h : a+i*h;
    yy = dSi(xx);
    ys1(1,i) = my_rectangles(xx,yy);
    ys2(1,i) = trapz(xx,yy);
    ys3(1,i) = my_simpson(xx,yy);
    %disp(i);
end
clear i xx yy;

figure(2);
plot(x,ys1, x,ys2, x,ys3);
legend({'my rectangles','trapz', 'my simpson'},'Location','southwest')
xlabel('x');
ylabel('y');
title('Integral of sin(x)/x');


%{
function y = dSi(x)
    y = sin(x)./x;
    y(1, find(x==0) ) = 1;
end
function I = my_rectangles(x,y)
    
    n = sum( size(x) ) - 2;
    a = x(1,1);
    b = x(1,end);
    h = (b-a)/n;
    
    
    I = 0;
    
    x1 = x(1,2:end);
    x2 = x(1,1:end-1);
    y1 = y(1,2:end);
    I = sum( y1(:).*(x1(:) - x2(:)) );
end
function I = my_simpson(x,y)
    
    n = sum( size(x) ) - 2;
    a = x(1,1);
    b = x(1,end);
    h = (b-a)/n;

    if mod(n+1,2) == 0
        %error('n is odd!');
        n = n+1;
        h = (b-a)/n;
        
        x = a:h:b; %n+1 elements
        y = dSi(x);%n+1 elements
        
        
        y1 = y(1,1:2:end-2);
        y2 = y(1,2:2:end-1);
        y3 = y(1,3:2:end);

        I = sum( (y1(:) + 4*y2(:) + y3(:)) );

        I = I*(h/3);
        return;
    end
    
    
    
    I = 0;
    
    y1 = y(1,1:2:end-2);
    y2 = y(1,2:2:end-1);
    y3 = y(1,3:2:end);
    
    I = sum( (y1(:) + 4*y2(:) + y3(:)) );
    
    I = I*(h/3);
    
end
%% task 12 (wait for check)
clear

a = -30;
b = 30;
n = 1000;

% Plot |h - h/2|
for i = 1:n
    
    h(i) = (b-a)/(2*i);
    
    xx = a : h(i) : b;
    yy = dSi(xx);
    
    
    
    Ihr = my_rectangles(xx,yy);
    Iht = trapz(xx,yy);
    Ihs = my_simpson(xx,yy);
    
    h2 = (b-a)/(4*i);
    xx2 = a : h2 : b;
    yy2 = dSi(xx2);
    
    Ih2r = my_rectangles(xx2,yy2);
    Ih2t = trapz(xx2,yy2);
    Ih2s = my_simpson(xx2,yy2);
    
    Ir(i) = abs(Ihr - Ih2r);
    It(i) = abs(Iht - Ih2t);
    Is(i) = abs(Ihs - Ih2s);
    
end

clear xx yy xx2 yy2 Ihr Iht Ihs Ih2r Ih2t Ih2s h2


figure(1);
loglog(h, Ir, h, It, h, Is);
title('Time difference');
legend({'|Ir-Ir2|','|It - It2|', '|Is - Is2|'},'Location','southwest')
xlabel('h');
ylabel('|dy - dy_n|');



%plot Si(x) by all methods

h = (b-a)/n;
x = a:h:b; %n+1 elements


ys1 = zeros(1,n+1);
ys2 = zeros(1,n+1);
ys3 = zeros(1,n+1);

for i = 1:n
    xx = a: h : a+i*h;
    yy = dSi(xx);
    ys1(1,i) = my_rectangles(xx,yy);
    ys2(1,i) = trapz(xx,yy);
    ys3(1,i) = my_simpson(xx,yy);
    %disp(i);
end
clear i xx yy;

figure(2);
plot(x,ys1, x,ys2, x,ys3);
legend({'my rectangles','trapz', 'my simpson'},'Location','southwest')
xlabel('x');
ylabel('y');
title('Integral of sin(x)/x');

%}


%% task 13 (complete)
clear


n = 500000;

h = logspace(-10,-1,n);


dot = 1;

dy = cos(dot);


dy_r = r_derivative(dot, h );
dy_c = c_derivative(dot, h );


diff_r_dotx(:) = abs( dy - dy_r( : ) );
diff_c_dotx(:) = abs( dy - dy_c( : ) );

loglog(h,diff_r_dotx, h, diff_c_dotx);
legend({'Right derivative','Central derivative'},'Location','southwest')
xlabel('h');
ylabel('|dy - dy_n|');


%{


%i = 2:n+1;
%h(:) = (b-a)./i(:);
diff_r_dotx = zeros(1,n);
diff_c_dotx = zeros(1,n);

diff_r_dotx(1,:) = abs( dy( (:) ) - dy_r( (:) ) );
diff_c_dotx(1,:) = abs( dy( (:) ) - dy_c( (:) ) );

diff_r_norm = zeros(1,n);
diff_c_norm = zeros(1,n);


h = zeros(1,n);
for i = 1:n
    
    h(i) = (b-a)/i;
    x = a:h(i):b;

    dy = cos(x);
    
    dy_r = r_derivative(x,h(i));
    dy_c = c_derivative(x,h(i));
    
    diff_r(i) = abs( dy(1)-dy_r(1) );
    diff_c(i) = abs( dy(1)-dy_r(1) );
    
    diff_r_norm(i) = norm(dy-dy_r);
    diff_c_norm(i) = norm(dy-dy_c);
end

loglog(h,diff_r_norm, h, diff_c_norm);

%}
%% training field
clear


error('oooooops.');

%% OLD 12
clc
clear
n = 3;
m = 4;
A = rand(n,m)

% Find G
i = 1:n;
j = 1:m/2;

G1(i(:),j(:)) = A(i(:),2*j(:))

%% ANIMATION
n = 50;
XY = 10 * rand(2,n) - 5;
for i=1:n
    plot(XY(1,i),XY(2,i),'or','MarkerSize',5,'MarkerFaceColor','r')
    axis([-5 5 -5 5])
    pause(.1)
end

