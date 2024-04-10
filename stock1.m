clear;
clc;

% ����������� ������ � ����� ��������� �������� 
p = 0.33;

% ���������� ��������� �������� (��� ��������� ����������� �������������� �. �.)
n = 20; 

% ����� �������� ����������� �������������� �. �.
count = 50000; 

% ��������� ������� ����������� �������������� ��������� ������� 
binomExample = binomrand(p, n, count);

% �������������� ������ � ������������ �������������� 
binDist = zeros(n + 1, 1);
for i=1:n+1
    binDist(i) = factorial(n) / (factorial(i-1)*factorial(n-i+1)) * p^(i-1) * (1-p)^(n-i+1);
end

% ������ ������ �������������
figure(1)
hold on;
histogram(binomExample, 'Normalization', 'probability');
bar(0:n, binDist, 'FaceColor', 'Green', 'BarWidth', 0.5);
legend('������������ �����������', '������������� �����������');
xlabel('\xi');
ylabel('P');
hold off;



% ���������� ������� ��������� ������� � �������������� ��������������
geomExample = geomrand(p, count);

% �������������� ������ � �������������� ��������������  
geomDist = zeros(max(geomExample) + 1, 1);
for i=1:max(geomExample)+1
    geomDist(i) = (1-p)^(i-1)*p;
end

% ������ �������������� �������������
figure(2)
hold on;
histogram(geomExample, 'Normalization', 'probability');
bar(0:max(geomExample), geomDist, 'FaceColor', 'Red', 'BarWidth', 0.45);
legend('������������ �����������', '������������� �����������');
xlabel('\xi');
ylabel('P');
hold off;

% �������� ���������� ������: P(X > t + s | X >= t) = P(X > s)
t = 2;

newGeomExample = sort(geomExample.*(geomExample>=t));
newGeomExample = newGeomExample(count - sum(geomExample>=t) + 1 : end);

figure(3)
hold on;
histogram(geomExample, 'Normalization', 'probability');
histogram(newGeomExample - t, 'Normalization', 'probability', 'FaceColor', 'Green', 'BinWidth', 0.45);
xlabel('\xi');
ylabel('P');
hold off;


% ���� � �������
N = 1000;
p = 0.5;
S = zeros(1, N);
Y= zeros(1, N);
for i=1:N
    if bernrand(p) == 1
        S(i) = S(i) + 1;
    else
        S(i) = S(i) - 1;
    end
    if i > 1
        S(i) = S(i) + S(i - 1);
    end
    Y(i) = S(i) / sqrt(N);
end
figure(4)
plot(1:N, Y)
title('���� � �������')
ylabel('Y(i)');
xlabel('i');



% ��������� ������� �. �. � �������������� ��������������
function x = geomrand(p, num)
    x = bernrand(p, 1, num);
    counter = 0;
    while nnz(x) < num
        x = x + (bernrand(p, 1, num) + x > zeros(1, num));
        counter = counter + 1;
    end
    x = ones(1, num) .* (counter + 1) - x;
end

% ��������� ������� �. �. � ������������ �������������� 
function x = binomrand(p, n, num)
    test_series = bernrand(p, n, num);
    x = sum(test_series);
end

% ��������� ����� ��������
function x = bernrand(p, varargin)
    if (p>1) || (p<0)   
        error('p = %4.2f is out of bounds [0, 1]',p)
    end
    x = rand(varargin{:}) < p;
end