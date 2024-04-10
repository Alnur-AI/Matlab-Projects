%% ������ 1 � 2
clear;
clc;

% ����� ������� ��������� �������� X (����� � ����� X - ��������� ��������� ��������)
n = 100;
% ����� ������� ��������� �������� Y = 1 - X
m = 90; 

alpha = 0.1; % ������� ����������
times = 1000; % ���������� ���������

counter_cant = 0; % ������� �������� �������� � ���, ��� ��������� �������� X �������� ����������
counter_symmetry = 0; % ������� �������� �������� � ���, ��� ��������� �������� X � Y=1-X ��������� ������������
counter_self_similarity = 0; % ������� �������� �������� � ���, ��������� �������� X/3 ������������ ��� ��, ��� � �. �. X ��� ������� X<=1/3


for i = 1:times
    % ���������� ������� � ������ �������� ������� ������������� ��������� �������� X
    [x, F] = cantrnd(n);
    % ��������� ���������� ����������� 
    empirF = sum((ones(n, n) .* x' < x), 2)./n;
    Dn = max(abs(empirF - F));
    % ��������� ���� ��������� ��������
    p = 1 - F_K(sqrt(n) * Dn);
    counter_cant = counter_cant + (p > alpha);
    
    % ���������� ������� ��������� �������� Y = 1 - X
    y = cantrnd(m);
    y = ones(m, 1) - y;
    % ��������� ���������� ��������
    empF_minus_empG = (sum((ones(n + m, n).* x' < cat(1,x,y)), 2) ./ n) - (sum((ones(n + m, m).* y' < cat(1,x,y)), 2) ./ m);
    Dnm = max(abs(empF_minus_empG));
    % ��������� ���� ��������� ��������
    p = 1 - F_K(sqrt(n*m/(n + m)) * Dnm);
    counter_symmetry = counter_symmetry + (p > alpha);
    
    % ���������� ������� ��������� ������� X ��� ������� X <= 1/3 � X/3
    m = sum(x <= 1/3); % ����� ������� ��������� �������� X ��� ������� X <= 1/3
    y = sort(x.*(x <= 1/3));
    y = y(n - m + 1 : end);
    x = x / 3;
    % ��������� ���������� ��������
    empF_minus_empG = (sum((ones(n + m, n).* x' < cat(1,x,y)), 2) ./ n) - (sum((ones(n + m, m).* y' < cat(1,x,y)), 2) ./ m);
    Dnm = max(abs(empF_minus_empG));
    % ��������� ���� ��������� ��������
    p = 1 - F_K(sqrt(n*m/(n + m)) * Dnm);
    counter_self_similarity = counter_self_similarity + (p > alpha);
end

str=sprintf('������� �������� �������� � ���, ��� ������� ������������� ��� ������ ������������ ��������� �������������, �����:');
disp(str);
counter_cant/times   

str=sprintf('������� �������� �������� � ���, ��� ��������� �������� X � Y = 1 - X ��������� ������������ (����� X - ��������� �. �.):');
disp(str);
counter_symmetry/times

str=sprintf('������� �������� �������� � ���, ��� ��������� ��������� �������� X �������� ��������� ����������� ������������ ������� �� 3:');
disp(str);
counter_self_similarity/times

%% �������� � ������� 1, 2

% ������������ �������� �����������
num = 10000;
[x, Fx] = cantrnd(num);
array_x = sortrows(cat(2, x, Fx));

arg_x = array_x(:,1);
expectedF_x = array_x(:,2);
empiricalF_x = sum(ones(num, num) .* arg_x' < arg_x, 2) / num;

figure;
hold on;
plot(arg_x, empiricalF_x)
plot(arg_x, expectedF_x)
legend('������������ ������� �������������', '������������� ������� �������������');
hold off; 

% ������������ �������� �������������� ������������ 0.5
y = cantrnd(num);
y = ones(num, 1) - y;
y = sort(y);

empiricalF_y = sum(ones(num, num) .* y' < y, 2) / num;

figure;
hold on;
plot(arg_x, empiricalF_x)
plot(y, empiricalF_y)
legend('������������ ������� ������������� �. �. X', '������������ ������� ������������� �. �. 1-X');
hold off;

% ������������ �������� �����������
x = cantrnd(num);
m = sum(x <= 1/3);
y = sort(x.*(x <= 1/3));
y = y(num - m + 1 : end);
x = sort(x) / 3;

empiricalF_y = sum(ones(m, m) .* y' < y, 2) / m;
empiricalF_x = sum(ones(num, num) .* x' < x, 2) / num;

figure;
hold on;
plot(x, empiricalF_x)
plot(y, empiricalF_y)
legend('������������ ������� ������������� �. �. X/3', '������������ ������� ������������� X ��� ������� X\leq1/3');
hold off;
%% ����� 3

N =1000;

sample_mean = zeros(N, 1);
sample_dispersion = zeros(N, 1);
for i = 1:N
    x = cantrnd(i);
    sample_mean(i) = (1/i) * sum(x);
    sample_dispersion(i) = (1/i) * sum( (x - ones(i,1).*sample_mean(i)).^2 );
end

figure(1);
hold on;
loglog(1:N, sample_mean);
plot(1:N, ones(N,1)*0.5);
legend('���������� �������', '�������������� ��������');
hold off;

figure(2);
hold on;
plot(1:N, sample_dispersion);
plot(1:N, ones(N,1)*0.125);
legend('���������� ���������', '���������');
hold off;

%% �������
% ������� ������������� �����������
function result = F_K(x, n)
    if nargin < 2
        n = 1000;
    end
    result = 1;
    for k = 1:n
        result = result + 2*((-1)^k)*exp(-2*(k^2)*(x^2));
    end
end

% ������� ��������� ���������� ��������� �������
function [x, F] = cantrnd(n, eps)
    if nargin < 1
        n = 1; % ����� �������
    end
    if nargin < 2
        eps = 1e-10; % �������� �������
    end
    m = round(-log(eps)/log(3));
    bern = bernrand(0.5, n, m);
    deg = -(1:m)';
    x = 2 * bern * 3.^deg;
    F = bern * 2.^deg;
end

% ������� ��������� ��������� ������� ��������
function x = bernrand(p, varargin)
    if (p>1) || (p<0)
        error('p = %4.2f is out of bounds [0, 1]',p)
    end
    x = rand(varargin{:}) < p;
end