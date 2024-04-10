%% Задание 6 основной пример
N = 40;
M = 30;
fHandle1 = @(x,y)(2-x.^3).*cos(x)-3*y.*exp(-y) + 2*cos(2*y);

mu1 = 1;
val = uNumerical(1,1,mu1,N,M);
x_ = linspace(0,1-1/N,N);
y_ = linspace(0,1-1/M,M);
[X,Y] =  meshgrid(y_,x_);
val_analit = uAnalytical(X,Y,1,1,mu1);
hold on
surf(X,Y,(real(val)),'FaceAlpha',0.5,'FaceColor','b','EdgeColor','none');
surf(X,Y,val_analit,'FaceAlpha',0.5,'FaceColor','r','EdgeColor','none');
%surf(X,Y,abs(val_analit-real(val)),'FaceAlpha',0.5,'EdgeColor','none');
xlabel('x');
ylabel('y');
zlabel('z');

%legend("Модуль разности");
legend("Численное решение","Аналитическое решение");

hold off 
view(3);
%% Задание 6 пример, когда ответ 3*sin(pi*x)+sin(pi*y)
N = 30;
M = 35;
mu1 = 1;
fHandle2 = @(x,y) -(pi^2+mu1)*(3*sin(pi*x)+sin(pi*y));
xiHandle2 = @(x) 3*sin(pi*x);
etaHandle2 = @(y) sin(pi*y);
% решение 3*sin(pix)+sin(piy);
val = solveDirichlet(fHandle2,xiHandle2,etaHandle2,mu1,N,M);
x_ = linspace(0,1-1/N,N);
y_ = linspace(0,1-1/M,M);
[X,Y] =  meshgrid(y_,x_);
f_1 =@(x,y) 3*sin(pi*x)+sin(y*pi);
val_f1 = f_1(X,Y);
hold on
%surf(X,Y,real(val),'FaceAlpha',0.5,'FaceColor','b');
%surf(X,Y,val_f1,'FaceAlpha',0.5,'FaceColor','r');
%legend('My solve','Real solve');
surf(X,Y,abs(real(val_f1)-(val)));
xlabel('x');
ylabel('y');
zlabel('z');
disp('Максимальная погрешность');
disp(max(max(abs(real(val)-val_f1))))
legend('abs (My solve-Real solve)');
hold off
view(3);
%% Задание 6 Пример, когда решение x(x-1)+y(y^3-1)
N = 30;
M = 40;
mu3 = 3;
fHandle3 = @(x,y) 2+12*y.^2-mu3*(x.*(x-1)+y.*(y.^3-1));
xiHandle3 = @(x) x.*(x-1);
etaHandle3 = @(y) y.*(y^3-1);
val3 = solveDirichlet(fHandle3,xiHandle3,etaHandle3,mu3,N,M);
x_ = linspace(0,1-1/N,N);
y_ = linspace(0,1-1/M,M);
[X,Y] =  meshgrid(y_,x_);
f_1 =@(x,y) x.*(x-1)+y.*(y.^3-1);
val_f1 = f_1(X,Y);
hold on
surf(X,Y,real(val3),'FaceAlpha',0.5,'FaceColor','b');
surf(X,Y,val_f1,'FaceAlpha',0.5,'FaceColor','r');
legend('My solve','Real solve');
%surf(X,Y,abs(real(val_f1)-(val3)));
%legend('abs (My solve-Real solve)');
disp('Максимальная погрешность');
disp(max(max(abs(real(val3)-val_f1))))
view(3);
xlabel('x');
ylabel('y');
zlabel('z');
hold off
%% Задание 6 примеры для отчета
N = 30;
M = 40;
mu4 = 3;
fHandle4 = @(x,y) exp(sin(x.*y))+sin(x.*exp(y));
xiHandle4 = @(x) sin(pi*x);
etaHandle4 = @(y) 4*y.*(y.^3-1) ;
x_ = linspace(0,1-1/N,N);
y_ = linspace(0,1-1/M,M);
[X,Y] =  meshgrid(y_,x_);
val4 = solveDirichlet(fHandle4,xiHandle4,etaHandle4,mu4,N,M);
surf(X,Y,real(val4),'FaceAlpha',0.5);
xlabel('x');
ylabel('y');
zlabel('z');
function res =  solveDirichlet(fHandle,xiHandle,etaHandle,mu,N,M)
    res = zeros(N,M);
    delta_x = 1/N;
    delta_y = 1/M;
    x_ = linspace(0,1-delta_x,N);
    y_ = linspace(0,1-delta_y,M);
    [X,Y] =  meshgrid(y_,x_);
        
    f = zeros(N,M);
    f = fHandle(X,Y);

    for (k=1:N)
        f(k,1)=0;
    end
    for (l=1:M)
        f(1,l)=0;
    end

    xi = zeros(1,M);
    for q = 1:M
        xi(q) = xiHandle(y_(q));
    end
    
    eta = zeros(1,N);
    for q = 1:N
        eta(q) = etaHandle(x_(q));
    end
    alpha = ifft(xi);
    beta = ifft(eta);
    C = zeros(N,M);
    for p=1:N
        for q=1:M
            C(p,q) = -4*(sin(pi*(p-1)/N)).^2/(delta_x.^2)-4*(sin(pi*(q-1)/M)).^2/(delta_y.^2)-mu;
        end
    end
    main_matrix = zeros(N+M-1,N+M-1);
    main_vector = zeros(1,N+M-1);
    %заполняем перввыми N числами
    D_mat = ifft2(f);
    for p=1:N
        % Вычислем значение правой части
        D = 0;
        for q=1:M
            D = D+(1/C(p,q))*D_mat(p,q);
        end
        main_vector(p) = beta(p)-D; 
        
        inv_c_vec = zeros(1,M);
        for i=1:M
            inv_c_vec(i)=1./C(p,i);
        end
        
        A = ifft(inv_c_vec)/N;
        % первый элемент заполняем отдельно
        for l=2:M
            main_matrix(p,l)=A(l);
        end
        % Заполняем первые M элементов p-ого столбца
        
        main_matrix(p,1) = main_matrix(p,1)+sum(inv_c_vec)./(M*N);
        % первый элемент дополняем отдельно
        
        for(k=2:N)
            main_matrix(p,k+M-1) = (sum(inv_c_vec)./(M*N))*exp((2*pi*1j*(k-1)*(p-1))/N);
        end
        %Вычисляем M+1 ...N+M-1 коэфф системы
        
        % Добавляем (p,1) последнее слагаемое
    end
    
    %Заполняем вторую часть матрицы
    for (q=2:M)
        D = 0;
        for p=1:N
            D = D+(1/C(p,q))*D_mat(p,q);
        end
        main_vector(q+N-1) = alpha(q)-D; 
        % вычисляем правую часть
        
        inv_c_vec = zeros(1,N);
        for i=1:N
            inv_c_vec(i)=1./C(i,q);
        end
        %вычислям вспомогательные значения
        
        main_matrix(q+N-1,1) = main_matrix(q+N-1,1)+sum(inv_c_vec)./(M*N);
        % Добавляем к 0(1) элменту
        
        for(l=2:M)
            main_matrix(q+N-1,l) = (sum(inv_c_vec)./(M*N))*exp((2*pi*1j*(l-1)*(q-1))/M);
        end
        % Заполняем первые M элементов
        
        B = ifft(inv_c_vec)/M;
        
        %Добавляем к 0(1) элементу
        for (k=2:N)
            main_matrix(q+N-1,M+k-1)=B(k);
        end
    end
    f_NM = bicg(main_matrix,main_vector',1e-7,1000);
    % Решаем систему
    for l=1:M
        f(1,l) = f_NM(l);
    end
    for k=2:N
        f(k,1) = f_NM(M+k-1);
    end % Заполняем матрицу f полученнымим значениями
    b_pq = ifft2(f);

    a_pq = zeros (N,M);
    a_pq = b_pq./C;
    res = real(fft2 (a_pq));
end
function res = uAnalytical(xMat,yMat,u1Zero,u2Zero,mu)
    val = zeros(size(xMat));
    u1 = @(s) exp(-(mu.^(1/2))*s)-exp((mu.^(1/2))*s);
    u2 = @(s) exp(2*mu.^(1/2)).*exp(-(mu.^(1/2))*s)-exp((mu.^(1/2))*s);
    t0 = 2*(mu.^(1/2))*(exp(2*mu.^(1/2))-1);
    
    f_x = @(s) (2-s.^3).*cos(s)+mu*u1Zero;
    f_y = @(s) -3.*s.*exp(-s)+cos(2*s)+mu*u2Zero;
    
    %f_x = @(s) 2*exp(2*s).*s.^2+0;
    %f_y = @(s) 0*s;
    N = size(xMat);
    M = size(xMat);
    for i=1:N(1)
        for j=1:M(2)
            x = xMat(i,j);
            y = yMat(i,j);
            set1 = linspace(0,x,1000);
            I1 = trapz(set1,u1(set1).*f_x(set1));
            I1 = u2(x).*I1/t0;
            
            set2 = linspace(x,1,1000);
            I2 = trapz(set2,u2(set2).*f_x(set2));
            I2 = u1(x).*I2/t0;
            
            set3 = linspace(0,y,1000);
            I3 = trapz(set3,u1(set3).*f_y(set3));
            I3 = u2(y).*I3/t0;
            
            set4 = linspace(y,1,1000);
            I4= trapz(set4,u2(set4).*f_y(set4));
            I4 = u1(y).*I4/t0;
            
            val(i,j) = I1+I2+I3+I4+u1Zero+u2Zero;
        end
    end
    res = val;
end
function res = uNumerical(u1Zero,u2Zero,mu_,N,M)
    fHangle = @(x,y) (2-x.^3).*cos(x)-3*y.*exp(-y)+cos(2*y);
    %fHangle = @(x,y) 2*x.^2.*exp(2*x)+0*y;
    
    xiHandle = @(x) uAnalytical(x,zeros(size(x)),u1Zero,u2Zero,mu_);
    etaHandle = @(y) uAnalytical(zeros(size(y)),y,u1Zero,u2Zero,mu_);
    res = solveDirichlet(fHangle,xiHandle,etaHandle,mu_,N,M);
end