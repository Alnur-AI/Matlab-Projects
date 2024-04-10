difference_and_errorEstimate(100);

function [] = difference_and_errorEstimate(N)
   big_constant = 10000;
   S = 0;
   for i = 1:big_constant
       S = S + ((-1)^i) / i^2;
       if i <= N
           Sn(i) = S;
       end
   end 
   grid = 1 : 1 : N;
   plot(grid, Sn - S);
   hold on;
   psi = ones(1, N) ./ grid;
   plot(grid, psi, 'r');
   hold off;
   legend ('S(n) - S', 'Оценка погрешности')
end