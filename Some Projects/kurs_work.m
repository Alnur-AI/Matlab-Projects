a = 0.7%2;
r = 6;%2;

Nstart = 0.01;%r*log(r)/(a*(r-1))+0.000016;
Pstart = 0.015;%log(r)/a;

iter_num = 25;
array = zeros(iter_num, 2);

array(1, :) = [Nstart, Pstart]
for i = 2:iter_num
    [array(i, 1), array(i, 2)] = f(a, r, array(i-1, 1), array(i-1, 2))
end

hold on;
plot(array(:, 1), array(:, 2), '-o', 'MarkerIndices', 1:iter_num);
%plot(array(2:7, 1), array(2:7, 2), 'o', 'MarkerIndices', 1, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'green', 'MarkerSize', 4);
plot(Nstart, Pstart, '-o', 'MarkerIndices', 1, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'red', 'MarkerSize', 4);
xlabel('N')
ylabel('P')
hold off;

%[a, b] = f(a, r, Nstart, Pstart)

function [Nnext, Pnext] = f(a, r, Ncurr, Pcurr)
    Nnext = r * Ncurr * exp(-a*Pcurr);
    Pnext = Ncurr * (1 - exp(-a*Pcurr));
end